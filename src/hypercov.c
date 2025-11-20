// HyperCov - Hyper Fast Coverage Calculation

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>
#include <sys/time.h>
#include <sam.h>
#include <hts.h>
#include <tbx.h>
#include <kstring.h>
#include <libdeflate.h>
#include <zlib.h>

// Add AVX2 support if available
#ifdef __AVX2__
#include <immintrin.h>
#define USE_AVX2 1
#else
#define USE_AVX2 0
#endif

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#define FORCE_INLINE __attribute__((always_inline)) inline
#define PREFETCH_READ(addr) __builtin_prefetch(addr, 0, 3)
#define PREFETCH_WRITE(addr) __builtin_prefetch(addr, 1, 3)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#define FORCE_INLINE inline
#define PREFETCH_READ(addr) ((void)0)
#define PREFETCH_WRITE(addr) ((void)0)
#endif
// ============================================================================
// CONSTANTS
// ============================================================================

#define WINDOW_SIZE         10000000
#define MIN_OVERLAP         150
#define MAX_OVERLAP         1000
#define BATCH_SIZE          16384
#define READ_BATCH_SIZE     256
#define COMPRESS_BUFFER_SIZE (4 << 20)
#define CACHE_LINE_SIZE     64
#define MAX_LINE_LENGTH     4096
#define SORT_BUFFER_SIZE    100000
#define CRAM_CACHE_SIZE     (1 << 20)
#define SINGLE_THREAD_BATCH_SIZE 512  // Larger batches for single thread

#define CACHE_ALIGNED __attribute__((aligned(CACHE_LINE_SIZE)))

// ============================================================================
// DATA STRUCTURES
// ============================================================================

typedef enum {
    MODE_WHOLE_GENOME,
    MODE_BED_REGIONS,
    MODE_GFF_GENES,
    MODE_FIXED_WINDOWS
} processing_mode_t;

typedef struct {
    char *input;
    char *output;
    char *reference;
    char *bed_file;
    char *gff_file;
    int threads;
    uint8_t min_mapq;
    uint16_t exclude_flags;
    uint32_t min_depth;
    uint64_t window_size;
    bool all_sites;
    bool unsorted;
    bool gc_content;
    processing_mode_t mode;
} args_t;

typedef struct {
    char *chrom;
    uint64_t start;
    uint64_t end;
    char *name;
    int tid;
} bed_region_t;

typedef struct {
    char *chrom;
    char *source;
    char *feature;
    uint64_t start;
    uint64_t end;
    char *strand;
    char *gene_id;
    char *gene_name;
    char *transcript_id;
    int tid;
} gff_feature_t;

typedef struct {
    int tid;
    char *chrom_name;
    uint64_t start;
    uint64_t end;
    char *name;
    uint64_t chrom_length;
    int region_idx;  // Direct index to avoid O(n) lookup
} region_t;

typedef struct {
    region_t *regions;
    size_t count;
    size_t capacity;
} region_list_t;

typedef struct {
    region_t *tasks;
    size_t total_tasks;
    CACHE_ALIGNED _Atomic size_t next_task;
} atomic_scheduler_t;

typedef struct {
    int32_t *depth_buffer;
    size_t buffer_size;
} simple_buffer_t;

typedef enum {
    MSG_BATCH_LINES,
    MSG_DONE
} writer_msg_type_t;

typedef struct {
    writer_msg_type_t type;
    char **lines;
    size_t count;
} writer_msg_t;

typedef struct {
    pthread_t thread;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    writer_msg_t *queue;
    size_t queue_size;
    size_t queue_capacity;
    size_t queue_head;
    size_t queue_tail;
    char *output_path;
    _Atomic bool running;
} writer_thread_t;

typedef struct {
    char *name;
    uint64_t total_bases;
    
    CACHE_ALIGNED _Atomic uint64_t covered_sites;
    CACHE_ALIGNED _Atomic int64_t total_depth;
    CACHE_ALIGNED _Atomic uint64_t gc_bases;
    CACHE_ALIGNED _Atomic uint64_t at_bases;
} region_stats_t;

typedef struct {
    bam1_t *record;
    int tid;
    int64_t pos;
} bam_record_sort_t;

typedef struct {
    int thread_id;
    samFile *fp;
    bam_hdr_t *hdr;
    hts_idx_t *idx;
    bam1_t **rec_batch;
    simple_buffer_t *buffer;
    atomic_scheduler_t *scheduler;
    args_t *args;
    region_stats_t *region_stats;
    size_t n_regions;
    writer_thread_t *perbase_writer;
    uint8_t **gc_data;
    bool is_cram;
    
    uint64_t *local_covered_sites;
    int64_t *local_total_depth;
    uint64_t *local_gc_bases;
    uint64_t *local_at_bases;
    
    CACHE_ALIGNED uint64_t local_reads_processed;
    CACHE_ALIGNED uint64_t local_reads_skipped;
} worker_context_t;

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

char* trim(char *str) {
    char *end;
    while (*str == ' ' || *str == '\t' || *str == '\n' || *str == '\r') str++;
    if (*str == 0) return str;
    end = str + strlen(str) - 1;
    while (end > str && (*end == ' ' || *end == '\t' || *end == '\n' || *end == '\r')) end--;
    *(end + 1) = 0;
    return str;
}

// ============================================================================
// BED FILE PARSING
// ============================================================================

region_list_t* parse_bed_file(const char *bed_file, bam_hdr_t *hdr) {
    region_list_t *list = (region_list_t*)calloc(1, sizeof(region_list_t));
    list->capacity = 10000;
    list->regions = (region_t*)malloc(list->capacity * sizeof(region_t));
    list->count = 0;
    
    gzFile fp = gzopen(bed_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open BED file: %s\n", bed_file);
        free(list->regions);
        free(list);
        return NULL;
    }
    
    char line[MAX_LINE_LENGTH];
    size_t line_num = 0;
    
    while (gzgets(fp, line, MAX_LINE_LENGTH)) {
        line_num++;
        
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        
        char chrom[256];
        uint64_t start, end;
        char name[256] = "";
        
        int n = sscanf(line, "%s\t%lu\t%lu\t%s", chrom, &start, &end, name);
        if (n < 3) {
            fprintf(stderr, "Warning: Skipping malformed BED line %zu\n", line_num);
            continue;
        }
        
        int tid = bam_name2id(hdr, chrom);
        if (tid < 0) {
            fprintf(stderr, "Warning: Chromosome '%s' not found in BAM header\n", chrom);
            continue;
        }
        
        if (list->count >= list->capacity) {
            list->capacity *= 2;
            list->regions = (region_t*)realloc(list->regions, 
                                              list->capacity * sizeof(region_t));
        }
        
        region_t *region = &list->regions[list->count];
        region->tid = tid;
        region->chrom_name = strdup(chrom);
        region->start = start;
        region->end = end;
        region->name = (name[0] != '\0') ? strdup(name) : NULL;
        region->chrom_length = hdr->target_len[tid];
        region->region_idx = list->count;  // Store index
        
        list->count++;
    }
    
    gzclose(fp);
    
    return list;
}

// ============================================================================
// GFF FILE PARSING
// ============================================================================

void parse_gff_attributes(char *attr_str, gff_feature_t *feature) {
    // GFF3 uses semicolons to separate, and = for key-value pairs
    // GTF uses semicolons and space-separated key "value" pairs
    
    char *saveptr = NULL;
    char *token = strtok_r(attr_str, ";", &saveptr);
    
    while (token) {
        token = trim(token);
        
        // Try GFF3 format first (key=value)
        char key[256], value[256];
        if (sscanf(token, "%[^=]=%s", key, value) == 2) {
            char *v = value;
            // Remove quotes if present
            if (*v == '"') v++;
            size_t len = strlen(v);
            if (len > 0 && v[len-1] == '"') v[len-1] = '\0';
            
            // Trim key as well
            char *k = trim(key);
            
            if (strcmp(k, "gene_id") == 0 || strcmp(k, "ID") == 0 || strcmp(k, "id") == 0) {
                feature->gene_id = strdup(v);
            } else if (strcmp(k, "gene_name") == 0 || strcmp(k, "Name") == 0 || strcmp(k, "name") == 0) {
                feature->gene_name = strdup(v);
            } else if (strcmp(k, "transcript_id") == 0 || strcmp(k, "Parent") == 0) {
                feature->transcript_id = strdup(v);
            }
        } 
        // Try GTF format (key "value")
        else if (sscanf(token, "%s \"%[^\"]\"", key, value) == 2) {
            char *k = trim(key);
            if (strcmp(k, "gene_id") == 0 || strcmp(k, "ID") == 0) {
                feature->gene_id = strdup(value);
            } else if (strcmp(k, "gene_name") == 0 || strcmp(k, "Name") == 0) {
                feature->gene_name = strdup(value);
            } else if (strcmp(k, "transcript_id") == 0) {
                feature->transcript_id = strdup(value);
            }
        }
        
        token = strtok_r(NULL, ";", &saveptr);
    }
}

region_list_t* parse_gff_file(const char *gff_file, bam_hdr_t *hdr, 
                               const char *feature_type) {
    region_list_t *list = (region_list_t*)calloc(1, sizeof(region_list_t));
    list->capacity = 50000;
    list->regions = (region_t*)malloc(list->capacity * sizeof(region_t));
    list->count = 0;
    
    gzFile fp = gzopen(gff_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open GFF file: %s\n", gff_file);
        free(list->regions);
        free(list);
        return NULL;
    }
    
    char line[MAX_LINE_LENGTH];
    size_t line_num = 0;
    size_t features_found = 0;
    size_t features_matched = 0;
    
    while (gzgets(fp, line, MAX_LINE_LENGTH)) {
        line_num++;
        
        // Skip comments and empty lines
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        
        char chrom[256], source[256], feature[256], strand[10];
        uint64_t start, end;
        char score[256], frame[256], attributes[2048];
        
        int n = sscanf(line, "%s\t%s\t%s\t%lu\t%lu\t%s\t%s\t%s\t%[^\n]",
                      chrom, source, feature, &start, &end, score, strand, frame, attributes);
        
        if (n < 9) continue;
        
        features_found++;
        
        // Case-insensitive comparison for feature type
        bool feature_matches = (strcasecmp(feature, feature_type) == 0);
        if (!feature_matches) continue;
        
        features_matched++;
        
        int tid = bam_name2id(hdr, chrom);
        if (tid < 0) {
            fprintf(stderr, "Warning: Chromosome '%s' not found in BAM header (line %zu)\n", chrom, line_num);
            continue;
        }
        
        gff_feature_t gff_feat = {0};
        gff_feat.chrom = strdup(chrom);
        gff_feat.source = strdup(source);
        gff_feat.feature = strdup(feature);
        gff_feat.start = start - 1;
        gff_feat.end = end;
        gff_feat.strand = strdup(strand);
        gff_feat.tid = tid;
        
        // Make a copy for parsing (strtok modifies the string)
        char attr_copy[2048];
        strncpy(attr_copy, attributes, sizeof(attr_copy) - 1);
        attr_copy[sizeof(attr_copy) - 1] = '\0';
        parse_gff_attributes(attr_copy, &gff_feat);
        
        if (list->count >= list->capacity) {
            list->capacity *= 2;
            list->regions = (region_t*)realloc(list->regions, 
                                              list->capacity * sizeof(region_t));
        }
        
        region_t *region = &list->regions[list->count];
        region->tid = tid;
        region->chrom_name = strdup(chrom);
        region->start = gff_feat.start;
        region->end = gff_feat.end;
        
        if (gff_feat.gene_name) {
            region->name = strdup(gff_feat.gene_name);
        } else if (gff_feat.gene_id) {
            region->name = strdup(gff_feat.gene_id);
        } else {
            char default_name[256];
            snprintf(default_name, sizeof(default_name), "%s:%lu-%lu", 
                    chrom, start, end);
            region->name = strdup(default_name);
        }
        
        region->chrom_length = hdr->target_len[tid];
        region->region_idx = list->count;  // Store index
        
        free(gff_feat.chrom);
        free(gff_feat.source);
        free(gff_feat.feature);
        free(gff_feat.strand);
        free(gff_feat.gene_id);
        free(gff_feat.gene_name);
        free(gff_feat.transcript_id);
        
        list->count++;
    }
    
    gzclose(fp);
    
    return list;
}

// ============================================================================
// WINDOW GENERATION
// ============================================================================

region_list_t* generate_processing_windows(bam_hdr_t *hdr) {
    region_list_t *list = (region_list_t*)calloc(1, sizeof(region_list_t));
    list->capacity = 10000;
    list->regions = (region_t*)malloc(list->capacity * sizeof(region_t));
    list->count = 0;
    
    for (int tid = 0; tid < hdr->n_targets; tid++) {
        uint64_t chrom_len = hdr->target_len[tid];
        
        for (uint64_t start = 0; start < chrom_len; start += WINDOW_SIZE) {
            uint64_t end = start + WINDOW_SIZE;
            if (end > chrom_len) end = chrom_len;
            
            if (list->count >= list->capacity) {
                list->capacity *= 2;
                list->regions = (region_t*)realloc(list->regions, 
                                                  list->capacity * sizeof(region_t));
            }
            
            region_t *region = &list->regions[list->count];
            region->tid = tid;
            region->chrom_name = strdup(hdr->target_name[tid]);
            region->start = start;
            region->end = end;
            
            char name[256];
            snprintf(name, sizeof(name), "%s:%lu-%lu", 
                    hdr->target_name[tid], start, end);
            region->name = strdup(name);
            region->chrom_length = chrom_len;
            region->region_idx = tid;  // Use tid for whole genome mode to aggregate per chromosome
            
            list->count++;
        }
    }
    
    return list;
}

region_list_t* generate_fixed_windows(bam_hdr_t *hdr, uint64_t window_size) {
    region_list_t *list = (region_list_t*)calloc(1, sizeof(region_list_t));
    list->capacity = 10000;
    list->regions = (region_t*)malloc(list->capacity * sizeof(region_t));
    list->count = 0;
    
    for (int tid = 0; tid < hdr->n_targets; tid++) {
        uint64_t chrom_len = hdr->target_len[tid];
        
        for (uint64_t start = 0; start < chrom_len; start += window_size) {
            uint64_t end = start + window_size;
            if (end > chrom_len) end = chrom_len;
            
            if (list->count >= list->capacity) {
                list->capacity *= 2;
                list->regions = (region_t*)realloc(list->regions, 
                                                  list->capacity * sizeof(region_t));
            }
            
            region_t *region = &list->regions[list->count];
            region->tid = tid;
            region->chrom_name = strdup(hdr->target_name[tid]);
            region->start = start;
            region->end = end;
            
            char name[256];
            snprintf(name, sizeof(name), "%s:%lu-%lu", 
                    hdr->target_name[tid], start, end);
            region->name = strdup(name);
            region->chrom_length = chrom_len;
            region->region_idx = list->count;
            
            list->count++;
        }
    }
    
    return list;
}

// Generate one region per chromosome for whole genome mode
region_list_t* generate_whole_chromosomes(bam_hdr_t *hdr) {
    region_list_t *list = (region_list_t*)calloc(1, sizeof(region_list_t));
    list->capacity = hdr->n_targets;
    list->regions = (region_t*)malloc(list->capacity * sizeof(region_t));
    list->count = 0;
    
    for (int tid = 0; tid < hdr->n_targets; tid++) {
        uint64_t chrom_len = hdr->target_len[tid];
        
        region_t *region = &list->regions[list->count];
        region->tid = tid;
        region->chrom_name = strdup(hdr->target_name[tid]);
        region->start = 0;
        region->end = chrom_len;
        region->name = strdup(hdr->target_name[tid]);
        region->chrom_length = chrom_len;
        region->region_idx = list->count;
        
        list->count++;
    }
    
    return list;
}

region_list_t* generate_genome_windows(bam_hdr_t *hdr) {
    return generate_fixed_windows(hdr, WINDOW_SIZE);
}

void free_region_list(region_list_t *list) {
    if (!list) return;
    
    for (size_t i = 0; i < list->count; i++) {
        free(list->regions[i].chrom_name);
        free(list->regions[i].name);
    }
    free(list->regions);
    free(list);
}

// ============================================================================
// ATOMIC SCHEDULER
// ============================================================================

atomic_scheduler_t* atomic_scheduler_create(region_t *tasks, size_t count) {
    atomic_scheduler_t *sched = (atomic_scheduler_t*)malloc(sizeof(atomic_scheduler_t));
    if (!sched) return NULL;
    
    sched->tasks = tasks;
    sched->total_tasks = count;
    atomic_init(&sched->next_task, 0);
    return sched;
}

static FORCE_INLINE bool atomic_scheduler_get_task(atomic_scheduler_t *sched, region_t *task) {
    size_t idx = atomic_fetch_add_explicit(&sched->next_task, 1, memory_order_relaxed);
    if (idx >= sched->total_tasks) return false;
    *task = sched->tasks[idx];
    return true;
}

void atomic_scheduler_free(atomic_scheduler_t *sched) {
    if (sched) free(sched);
}

// ============================================================================
// SIMPLE BUFFER
// ============================================================================

simple_buffer_t* simple_buffer_create(size_t size) {
    simple_buffer_t *buf = (simple_buffer_t*)malloc(sizeof(simple_buffer_t));
    if (!buf) return NULL;
    
    buf->buffer_size = size;
    // Use cache-aligned allocation
    buf->depth_buffer = (int32_t*)aligned_alloc(CACHE_LINE_SIZE, size * sizeof(int32_t));
    if (!buf->depth_buffer) {
        free(buf);
        return NULL;
    }
    memset(buf->depth_buffer, 0, size * sizeof(int32_t));
    return buf;
}

void simple_buffer_free(simple_buffer_t *buf) {
    if (!buf) return;
    free(buf->depth_buffer);  // aligned_alloc uses regular free()
    free(buf);
}

// ============================================================================
// WRITER THREAD
// ============================================================================

void* writer_thread_func(void *arg) {
    writer_thread_t *writer = (writer_thread_t*)arg;
    
    FILE *fp = fopen(writer->output_path, "wb");
    if (!fp) return NULL;
    
    struct libdeflate_compressor *comp = libdeflate_alloc_compressor(6);
    if (!comp) {
        fclose(fp);
        return NULL;
    }
    
    char *buffer = (char*)malloc(COMPRESS_BUFFER_SIZE);
    char *compressed = (char*)malloc(COMPRESS_BUFFER_SIZE);
    if (!buffer || !compressed) {
        free(buffer);
        free(compressed);
        libdeflate_free_compressor(comp);
        fclose(fp);
        return NULL;
    }
    
    size_t buffer_used = 0;
    
    while (atomic_load_explicit(&writer->running, memory_order_acquire)) {
        pthread_mutex_lock(&writer->mutex);
        
        while (writer->queue_size == 0 && 
               atomic_load_explicit(&writer->running, memory_order_acquire)) {
            pthread_cond_wait(&writer->cond, &writer->mutex);
        }
        
        if (writer->queue_size == 0 && 
            !atomic_load_explicit(&writer->running, memory_order_acquire)) {
            pthread_mutex_unlock(&writer->mutex);
            break;
        }
        
        writer_msg_t msg = writer->queue[writer->queue_head];
        writer->queue_head = (writer->queue_head + 1) % writer->queue_capacity;
        writer->queue_size--;
        
        pthread_mutex_unlock(&writer->mutex);
        
        if (msg.type == MSG_DONE) break;
        
        for (size_t i = 0; i < msg.count; i++) {
            size_t len = strlen(msg.lines[i]);
            
            if (buffer_used + len + 1 >= COMPRESS_BUFFER_SIZE - 4096) {
                size_t comp_size = libdeflate_gzip_compress(comp, buffer, buffer_used,
                                                           compressed, COMPRESS_BUFFER_SIZE);
                if (comp_size > 0) fwrite(compressed, 1, comp_size, fp);
                buffer_used = 0;
            }
            
            memcpy(buffer + buffer_used, msg.lines[i], len);
            buffer_used += len;
            buffer[buffer_used++] = '\n';
            free(msg.lines[i]);
        }
        
        free(msg.lines);
    }
    
    if (buffer_used > 0) {
        size_t comp_size = libdeflate_gzip_compress(comp, buffer, buffer_used,
                                                    compressed, COMPRESS_BUFFER_SIZE);
        if (comp_size > 0) fwrite(compressed, 1, comp_size, fp);
    }
    
    free(buffer);
    free(compressed);
    libdeflate_free_compressor(comp);
    fclose(fp);
    return NULL;
}

writer_thread_t* writer_thread_create(const char *output_path) {
    writer_thread_t *writer = (writer_thread_t*)calloc(1, sizeof(writer_thread_t));
    if (!writer) return NULL;
    
    writer->output_path = strdup(output_path);
    if (!writer->output_path) {
        free(writer);
        return NULL;
    }
    
    writer->queue_capacity = 1024;
    writer->queue = (writer_msg_t*)calloc(writer->queue_capacity, sizeof(writer_msg_t));
    if (!writer->queue) {
        free(writer->output_path);
        free(writer);
        return NULL;
    }
    
    writer->queue_size = 0;
    writer->queue_head = 0;
    writer->queue_tail = 0;
    atomic_init(&writer->running, true);
    pthread_mutex_init(&writer->mutex, NULL);
    pthread_cond_init(&writer->cond, NULL);
    
    if (pthread_create(&writer->thread, NULL, writer_thread_func, writer) != 0) {
        free(writer->queue);
        free(writer->output_path);
        pthread_mutex_destroy(&writer->mutex);
        pthread_cond_destroy(&writer->cond);
        free(writer);
        return NULL;
    }
    
    return writer;
}

void writer_thread_send(writer_thread_t *writer, writer_msg_t msg) {
    if (!writer) {
        if (msg.type == MSG_BATCH_LINES && msg.lines) {
            for (size_t i = 0; i < msg.count; i++) {
                free(msg.lines[i]);
            }
            free(msg.lines);
        }
        return;
    }
    
    pthread_mutex_lock(&writer->mutex);
    
    while (writer->queue_size >= writer->queue_capacity - 1) {
        pthread_mutex_unlock(&writer->mutex);
        usleep(50);
        pthread_mutex_lock(&writer->mutex);
    }
    
    writer->queue[writer->queue_tail] = msg;
    writer->queue_tail = (writer->queue_tail + 1) % writer->queue_capacity;
    writer->queue_size++;
    
    pthread_cond_signal(&writer->cond);
    pthread_mutex_unlock(&writer->mutex);
}

void writer_thread_close(writer_thread_t *writer) {
    if (!writer) return;
    
    // Ensure all queued messages are processed before sending DONE
    pthread_mutex_lock(&writer->mutex);
    while (writer->queue_size > 0) {
        pthread_mutex_unlock(&writer->mutex);
        usleep(1000);  // Wait for queue to drain
        pthread_mutex_lock(&writer->mutex);
    }
    pthread_mutex_unlock(&writer->mutex);
    
    writer_msg_t msg = { .type = MSG_DONE, .lines = NULL, .count = 0 };
    writer_thread_send(writer, msg);
    
    atomic_store_explicit(&writer->running, false, memory_order_release);
    pthread_cond_signal(&writer->cond);
    pthread_join(writer->thread, NULL);
    
    for (size_t i = 0; i < writer->queue_size; i++) {
        size_t idx = (writer->queue_head + i) % writer->queue_capacity;
        if (writer->queue[idx].type == MSG_BATCH_LINES && writer->queue[idx].lines) {
            for (size_t j = 0; j < writer->queue[idx].count; j++) {
                free(writer->queue[idx].lines[j]);
            }
            free(writer->queue[idx].lines);
        }
    }
    
    pthread_mutex_destroy(&writer->mutex);
    pthread_cond_destroy(&writer->cond);
    free(writer->output_path);
    free(writer->queue);
    free(writer);
}

// ============================================================================
// UNSORTED BAM HANDLING
// ============================================================================

int compare_bam_records(const void *a, const void *b) {
    const bam_record_sort_t *ra = (const bam_record_sort_t*)a;
    const bam_record_sort_t *rb = (const bam_record_sort_t*)b;
    
    if (ra->tid != rb->tid) return ra->tid - rb->tid;
    if (ra->pos != rb->pos) return (ra->pos > rb->pos) ? 1 : -1;
    return 0;
}

char* sort_bam_file(const char *input, const char *reference, int threads) {
    (void)threads;  // Unused parameter, but kept for API compatibility
    fprintf(stderr, "Sorting unsorted BAM file (this may take a while)...\n");
    double start_time = get_time();
    
    char *temp_file = (char*)malloc(256);
    snprintf(temp_file, 256, "%s.sorted.tmp.bam", input);
    
    samFile *in = sam_open(input, "r");
    if (!in) {
        fprintf(stderr, "Error: Cannot open input BAM\n");
        free(temp_file);
        return NULL;
    }
    
    if (reference) hts_set_fai_filename(in, reference);
    
    bam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
        fprintf(stderr, "Error: Cannot read header\n");
        sam_close(in);
        free(temp_file);
        return NULL;
    }
    
    fprintf(stderr, "Loading BAM records into memory...\n");
    bam_record_sort_t *records = (bam_record_sort_t*)malloc(
        SORT_BUFFER_SIZE * sizeof(bam_record_sort_t));
    size_t n_records = 0;
    size_t capacity = SORT_BUFFER_SIZE;
    
    bam1_t *rec = bam_init1();
    while (sam_read1(in, hdr, rec) >= 0) {
        if (n_records >= capacity) {
            capacity *= 2;
            records = (bam_record_sort_t*)realloc(records, 
                capacity * sizeof(bam_record_sort_t));
        }
        
        records[n_records].record = bam_dup1(rec);
        records[n_records].tid = rec->core.tid;
        records[n_records].pos = rec->core.pos;
        n_records++;
        
        if (n_records % 1000000 == 0) {
            fprintf(stderr, "  Loaded %zu million records...\n", n_records / 1000000);
        }
    }
    
    bam_destroy1(rec);
    sam_close(in);
    
    fprintf(stderr, "Loaded %zu records, sorting...\n", n_records);
    qsort(records, n_records, sizeof(bam_record_sort_t), compare_bam_records);
    
    fprintf(stderr, "Writing sorted BAM...\n");
    
    samFile *out = sam_open(temp_file, "wb");
    if (!out) {
        fprintf(stderr, "Error: Cannot create temporary file\n");
        free(temp_file);
        return NULL;
    }
    
    sam_hdr_write(out, hdr);
    
    for (size_t i = 0; i < n_records; i++) {
        sam_write1(out, hdr, records[i].record);
        bam_destroy1(records[i].record);
        
        if ((i + 1) % 1000000 == 0) {
            fprintf(stderr, "  Wrote %zu million records...\n", (i + 1) / 1000000);
        }
    }
    
    sam_close(out);
    bam_hdr_destroy(hdr);
    free(records);
    
    fprintf(stderr, "Indexing sorted BAM...\n");
    sam_index_build(temp_file, 0);
    
    double elapsed = get_time() - start_time;
    fprintf(stderr, "Sorting completed in %.2f seconds\n", elapsed);
    
    return temp_file;
}

// ============================================================================
// PAIRED-END OVERLAP HANDLING
// ============================================================================

// Identify if this record is the leftmost mate of a proper pair on the same contig.
static inline bool is_leftmost_pair(const bam1_t *rec) {
    if (!(rec->core.flag & BAM_FPAIRED) ||
        !(rec->core.flag & BAM_FPROPER_PAIR) ||
        rec->core.tid != rec->core.mtid) {
        return false;
    }
    if (rec->core.pos < rec->core.mpos) return true;
    if (rec->core.pos > rec->core.mpos) return false;
    // Tie-break when starts equal: treat READ1 as leftmost
    return (rec->core.flag & BAM_FREAD1) != 0;
}

static inline bool should_skip_read(bam1_t *rec) {
    // Never skip the right mate; overlapping handling is done via truncation on the leftmost mate only.
    // Keeping this function for call-site compatibility.
    (void)rec;
    return false;
}

static inline int64_t get_truncate_position(bam1_t *rec) {
    // Truncate only the leftmost read at the mate start when there is an actual overlap.
    if (!is_leftmost_pair(rec)) return -1;

    int64_t read_end = bam_endpos(rec);
    int64_t mate_start = rec->core.mpos;

    if (read_end > mate_start) {
        return mate_start;  // Truncate leftmost read at mate start to avoid double counting
    }
    return -1;
}

// ============================================================================
// VECTORIZED COVERAGE OPERATIONS
// ============================================================================

#if USE_AVX2
static FORCE_INLINE void add_coverage_avx2(int32_t *depths,
                                           uint64_t idx_start,
                                           uint64_t idx_end) {
    uint64_t i = idx_start;
    __m256i ones = _mm256_set1_epi32(1);
    
    // Process 8 int32_t at a time
    for (; i + 8 <= idx_end; i += 8) {
        __m256i depth_vec = _mm256_loadu_si256((__m256i*)&depths[i]);
        depth_vec = _mm256_add_epi32(depth_vec, ones);
        _mm256_storeu_si256((__m256i*)&depths[i], depth_vec);
    }
    
    // Remainder
    for (; i < idx_end; i++) {
        depths[i]++;
    }
}

static inline uint64_t horizontal_sum_epi32(__m256i vec) {
    __m128i low = _mm256_castsi256_si128(vec);
    __m128i high = _mm256_extracti128_si256(vec, 1);
    __m128i sum128 = _mm_add_epi32(low, high);
    __m128i hi64 = _mm_unpackhi_epi64(sum128, sum128);
    __m128i sum64 = _mm_add_epi32(sum128, hi64);
    __m128i hi32 = _mm_shuffle_epi32(sum64, _MM_SHUFFLE(2, 3, 0, 1));
    __m128i sum32 = _mm_add_epi32(sum64, hi32);
    return _mm_cvtsi128_si32(sum32);
}

static void compute_stats_avx2(int32_t *depths, uint64_t len, int32_t min_depth,
                               uint64_t *covered, int64_t *total) {
    __m256i min_vec = _mm256_set1_epi32(min_depth - 1);
    __m256i cov_vec = _mm256_setzero_si256();
    __m256i tot_vec = _mm256_setzero_si256();
    
    uint64_t i = 0;
    for (; i + 8 <= len; i += 8) {
        __m256i depth_vec = _mm256_loadu_si256((__m256i*)&depths[i]);
        
        // Compare: depth >= min_depth
        __m256i cmp = _mm256_cmpgt_epi32(depth_vec, min_vec);
        __m256i ones = _mm256_srli_epi32(cmp, 31);
        cov_vec = _mm256_add_epi32(cov_vec, ones);
        
        // Sum depths
        tot_vec = _mm256_add_epi32(tot_vec, depth_vec);
    }
    
    // Horizontal sum
    *covered = horizontal_sum_epi32(cov_vec);
    *total = horizontal_sum_epi32(tot_vec);
    
    // Remainder
    for (; i < len; i++) {
        int32_t d = depths[i];
        *total += d;
        if (d >= min_depth) (*covered)++;
    }
}
#endif

// ============================================================================
// FAST CIGAR PROCESSING
// ============================================================================

// Pre-computed lookup tables for CIGAR operations
static const uint8_t cigar_is_match[16] = {
    [BAM_CMATCH] = 1, [BAM_CEQUAL] = 1, [BAM_CDIFF] = 1,
};

static const uint8_t cigar_consumes_ref[16] = {
    [BAM_CMATCH] = 1, [BAM_CEQUAL] = 1, [BAM_CDIFF] = 1,
    [BAM_CDEL] = 1, [BAM_CREF_SKIP] = 1,
};

// Fast path: single match operation (most common case)
static FORCE_INLINE void add_read_simple(int32_t read_start, int32_t read_end,
                                         int32_t *depths, uint64_t win_start, 
                                         uint64_t win_end) {
    // Clip to window
    if (read_start < (int32_t)win_start) read_start = win_start;
    if (read_end > (int32_t)win_end) read_end = win_end;
    
    if (read_start >= read_end) return;
    
    uint64_t idx_start = read_start - win_start;
    uint64_t idx_end = read_end - win_start;
    
#if USE_AVX2
    // Use AVX2 for long spans
    if (idx_end - idx_start > 32) {
        add_coverage_avx2(depths, idx_start, idx_end);
    } else {
        for (uint64_t j = idx_start; j < idx_end; j++) {
            depths[j]++;
        }
    }
#else
    for (uint64_t j = idx_start; j < idx_end; j++) {
        depths[j]++;
    }
#endif
}

// Medium path: unrolled loop for 2-4 operations
static FORCE_INLINE void process_cigar_unrolled(bam1_t *rec, const uint32_t *cigar,
                                                 uint32_t n_cigar, int32_t *depths,
                                                 uint64_t win_start, uint64_t win_end,
                                                 int64_t truncate_pos) {
    int32_t pos = rec->core.pos;
    
    // Manually unroll for common cases (2-4 operations)
    for (uint32_t i = 0; i < n_cigar; i++) {
        const int op = bam_cigar_op(cigar[i]);
        const int op_len = bam_cigar_oplen(cigar[i]);
        
        if (cigar_is_match[op]) {
            int32_t start = pos;
            int32_t end = pos + op_len;
            
            // Apply truncation
            if (truncate_pos > 0 && end > truncate_pos) {
                end = truncate_pos;
            }
            
            if (start < end) {
                add_read_simple(start, end, depths, win_start, win_end);
            }
            
            pos += op_len;
            
            if (truncate_pos > 0 && pos >= truncate_pos) break;
        } else if (cigar_consumes_ref[op]) {
            pos += op_len;
        }
    }
}

// Standard path: full CIGAR processing
static FORCE_INLINE void process_cigar_standard(bam1_t *rec, const uint32_t *cigar,
                                                 uint32_t n_cigar, int32_t *depths,
                                                 uint64_t win_start, uint64_t win_end,
                                                 int64_t truncate_pos) {
    int32_t pos = rec->core.pos;
    
    for (uint32_t i = 0; i < n_cigar; i++) {
        const int op = bam_cigar_op(cigar[i]);
        const int op_len = bam_cigar_oplen(cigar[i]);
        
        if (cigar_is_match[op]) {
            int32_t start = pos;
            int32_t end = pos + op_len;
            
            if (truncate_pos > 0 && end > truncate_pos) {
                end = truncate_pos;
            }
            
            if (start < end) {
                add_read_simple(start, end, depths, win_start, win_end);
            }
            
            pos += op_len;
            
            if (truncate_pos > 0 && pos >= truncate_pos) break;
        } else if (cigar_consumes_ref[op]) {
            pos += op_len;
        }
    }
}

// Main optimized alignment processor with fast paths
static FORCE_INLINE void process_alignment_optimized(bam1_t *rec, int32_t *depths,
                                                      uint64_t win_start, uint64_t win_end,
                                                      int64_t truncate_pos) {
    const uint32_t *cigar = bam_get_cigar(rec);
    const uint32_t n_cigar = rec->core.n_cigar;
    
    // FAST PATH: Single operation (e.g., "150M") - handles ~70% of reads
    if (LIKELY(n_cigar == 1)) {
        const uint32_t op = bam_cigar_op(cigar[0]);
        if (LIKELY(cigar_is_match[op])) {
            int32_t start = rec->core.pos;
            int32_t end = start + bam_cigar_oplen(cigar[0]);
            
            if (truncate_pos > 0 && end > truncate_pos) {
                end = truncate_pos;
            }
            
            add_read_simple(start, end, depths, win_start, win_end);
            return;
        }
    }
    
    // MEDIUM PATH: 2-4 operations - handles ~25% of reads
    if (LIKELY(n_cigar <= 4)) {
        process_cigar_unrolled(rec, cigar, n_cigar, depths, win_start, win_end, truncate_pos);
        return;
    }
    
    // SLOW PATH: Complex CIGAR - handles ~5% of reads
    process_cigar_standard(rec, cigar, n_cigar, depths, win_start, win_end, truncate_pos);
}

// ============================================================================
// GC CONTENT AND HELPER FUNCTIONS
// ============================================================================

uint8_t** load_gc_content(const char *fasta_path, bam_hdr_t *hdr) {
    gzFile fp = gzopen(fasta_path, "r");
    if (!fp) {
        fprintf(stderr, "Warning: Cannot open reference FASTA: %s\n", fasta_path);
        return NULL;
    }
    
    uint8_t **gc_data = (uint8_t**)calloc(hdr->n_targets, sizeof(uint8_t*));
    if (!gc_data) {
        gzclose(fp);
        return NULL;
    }
    
    char line[8192];
    int current_tid = -1;
    size_t current_pos = 0;
    uint8_t *current_seq = NULL;
    
    while (gzgets(fp, line, sizeof(line))) {
        if (line[0] == '>') {
            char chrom[256];
            if (sscanf(line + 1, "%s", chrom) == 1) {
                current_tid = bam_name2id(hdr, chrom);
                if (current_tid >= 0) {
                    current_seq = (uint8_t*)malloc(hdr->target_len[current_tid]);
                    if (current_seq) {
                        gc_data[current_tid] = current_seq;
                        current_pos = 0;
                    } else {
                        current_tid = -1;
                    }
                } else {
                    current_seq = NULL;
                }
            }
        } else if (current_seq && current_tid >= 0) {
            for (char *p = line; *p && *p != '\n' && *p != '\r'; p++) {
                if (current_pos >= hdr->target_len[current_tid]) break;
                uint8_t gc_flag;
                switch (*p) {
                    case 'G': case 'g': case 'C': case 'c':
                        gc_flag = 2; break;
                    case 'A': case 'a': case 'T': case 't':
                        gc_flag = 1; break;
                    default:
                        gc_flag = 0; break;
                }
                current_seq[current_pos++] = gc_flag;
            }
        }
    }
    
    gzclose(fp);
    return gc_data;
}

void free_gc_data(uint8_t **gc_data, int n_targets) {
    if (!gc_data) return;
    for (int i = 0; i < n_targets; i++) {
        free(gc_data[i]);
    }
    free(gc_data);
}

static FORCE_INLINE void calculate_gc_content(
    uint8_t *gc_seq,
    uint64_t start,
    uint64_t end,
    uint64_t *gc_bases,
    uint64_t *at_bases)
{
    uint64_t gc = 0, at = 0;
    for (uint64_t i = start; i < end; i++) {
        gc += (gc_seq[i] == 2);
        at += (gc_seq[i] == 1);
    }
    *gc_bases = gc;
    *at_bases = at;
}

static void configure_input_file(samFile *fp, bool is_cram, const char *reference) {
    if (is_cram) {
        hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0);
        hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, 
                   SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
        hts_set_opt(fp, HTS_OPT_CACHE_SIZE, CRAM_CACHE_SIZE);
        hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, 256 * 1024);
        if (reference) {
            hts_set_fai_filename(fp, reference);
        }
    } else {
        hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, 256 * 1024);
    }
    hts_set_threads(fp, 1);
    hts_set_opt(fp, HTS_OPT_THREAD_POOL, NULL);
}

static bool is_cram_file(const char *filename) {
    size_t len = strlen(filename);
    return (len > 5 && strcmp(filename + len - 5, ".cram") == 0) ||
           (len > 8 && strcmp(filename + len - 8, ".cram.gz") == 0);
}

// ============================================================================
// COVERAGE CALCULATION
// ============================================================================

int process_region(region_t *task, worker_context_t *ctx) {
    const uint64_t win_start = task->start;
    const uint64_t win_end = task->end;
    const uint64_t win_len = win_end - win_start;
    const int tid = task->tid;
    const int task_idx = task->region_idx;
    
    if (win_len > ctx->buffer->buffer_size) return -1;
    
    int32_t *depths = ctx->buffer->depth_buffer;
    memset(depths, 0, win_len * sizeof(int32_t));
    
    hts_itr_t *iter = sam_itr_queryi(ctx->idx, tid, win_start, win_end);
    if (!iter) return -1;
    
    const uint32_t flags = ctx->args->exclude_flags;
    const uint8_t min_mapq = ctx->args->min_mapq;
    bam1_t **rec_batch = ctx->rec_batch;
    
    uint64_t reads_processed = 0;
    uint64_t reads_skipped = 0;
    
    // Simple batched reading
    while (true) {
        int batch_count = 0;
        for (int b = 0; b < 128; b++) {
            if (sam_itr_next(ctx->fp, iter, rec_batch[b]) < 0) break;
            batch_count++;
        }
        if (batch_count == 0) break;
        
        // Process batch with optimized fast-path
        for (int b = 0; b < batch_count; b++) {
            bam1_t *rec = rec_batch[b];
            
            if ((rec->core.flag & flags) || (rec->core.flag & BAM_FUNMAP) || rec->core.qual < min_mapq) {
                reads_skipped++;
                continue;
            }
            
            // Skip overlapping paired-end reads
            if (should_skip_read(rec)) {
                reads_skipped++;
                continue;
            }
            
            // REMOVE TRUNCATION to match PanDepth behavior
            // int64_t truncate_pos = get_truncate_position(rec);
            int64_t truncate_pos = -1;  // Disable truncation
            
            // Use optimized fast-path processing
            process_alignment_optimized(rec, depths, win_start, win_end, truncate_pos);
            
            reads_processed++;
        }
    }
    
    hts_itr_destroy(iter);
    
    ctx->local_reads_processed += reads_processed;
    ctx->local_reads_skipped += reads_skipped;
    
    // Compute statistics with vectorization
    uint64_t covered = 0;
    int64_t total = 0;
    const int32_t min_depth = ctx->args->min_depth;
    
#if USE_AVX2
    compute_stats_avx2(depths, win_len, min_depth, &covered, &total);
#else
    for (uint64_t i = 0; i < win_len; i++) {
        int32_t d = depths[i];
        total += d;
        if (d >= min_depth) covered++;
    }
#endif
    
    ctx->local_covered_sites[task_idx] += covered;
    ctx->local_total_depth[task_idx] += total;
    
    if (ctx->args->gc_content && ctx->gc_data && ctx->gc_data[tid]) {
        uint64_t gc_bases = 0;
        uint64_t at_bases = 0;
        calculate_gc_content(ctx->gc_data[tid], win_start, win_end, &gc_bases, &at_bases);
        ctx->local_gc_bases[task_idx] += gc_bases;
        ctx->local_at_bases[task_idx] += at_bases;
    }
    
    return 0;
}

int process_region_single_thread(region_t *task, samFile *fp, bam_hdr_t *hdr, 
                                  hts_idx_t *idx, bam1_t **rec_batch, 
                                  int32_t *depths, args_t *args,
                                  uint64_t *covered_out, int64_t *total_out,
                                  uint64_t *gc_out, uint64_t *at_out,
                                  uint8_t **gc_data) {
    (void)hdr;  // Unused parameter - kept for API consistency
    const uint64_t win_start = task->start;
    const uint64_t win_end = task->end;
    const uint64_t win_len = win_end - win_start;
    const int tid = task->tid;
    
    // Zero buffer - use direct pointer for speed
    memset(depths, 0, win_len * sizeof(int32_t));
    
    hts_itr_t *iter = sam_itr_queryi(idx, tid, win_start, win_end);
    if (!iter) return -1;
    
    const uint32_t flags = args->exclude_flags;
    const uint8_t min_mapq = args->min_mapq;
    
    // Larger batches for single thread
    while (true) {
        int batch_count = 0;
        for (int b = 0; b < SINGLE_THREAD_BATCH_SIZE; b++) {
            if (sam_itr_next(fp, iter, rec_batch[b]) < 0) break;
            batch_count++;
        }
        if (batch_count == 0) break;
        
        // Tight loop - no function calls in hot path
        for (int b = 0; b < batch_count; b++) {
            bam1_t *rec = rec_batch[b];
            
            // Quick filter
            if ((rec->core.flag & (args->exclude_flags | BAM_FUNMAP)) || rec->core.qual < args->min_mapq) {
                continue;
            }
            
            // Use common helpers for paired-end overlap handling
            if (should_skip_read(rec)) {
                continue;
            }
            // REMOVE TRUNCATION to match PanDepth behavior
            // int64_t truncate_pos = get_truncate_position(rec);
            int64_t truncate_pos = -1;  // Disable truncation
            
            // Inline CIGAR processing for maximum speed
            const uint32_t *cigar = bam_get_cigar(rec);
            const uint32_t n_cigar = rec->core.n_cigar;
            
            // Fast path: single match operation
            if (LIKELY(n_cigar == 1 && cigar_is_match[bam_cigar_op(cigar[0])])) {
                int32_t start = rec->core.pos;
                int32_t end = start + bam_cigar_oplen(cigar[0]);
                
                if (truncate_pos > 0 && end > truncate_pos) end = truncate_pos;
                if (start < (int32_t)win_start) start = win_start;
                if (end > (int32_t)win_end) end = win_end;
                
                if (start < end) {
                    uint64_t idx_start = start - win_start;
                    uint64_t idx_end = end - win_start;
                    
                    // Direct increment - fastest path
                    for (uint64_t j = idx_start; j < idx_end; j++) {
                        depths[j]++;
                    }
                }
            } else {
                // Standard CIGAR processing
                int32_t pos = rec->core.pos;
                for (uint32_t i = 0; i < n_cigar; i++) {
                    const int op = bam_cigar_op(cigar[i]);
                    const int op_len = bam_cigar_oplen(cigar[i]);
                    
                    if (cigar_is_match[op]) {
                        int32_t start = pos;
                        int32_t end = pos + op_len;
                        
                        if (truncate_pos > 0 && end > truncate_pos) end = truncate_pos;
                        if (start < (int32_t)win_start) start = win_start;
                        if (end > (int32_t)win_end) end = win_end;
                        
                        if (start < end) {
                            uint64_t idx_start = start - win_start;
                            uint64_t idx_end = end - win_start;
                            for (uint64_t j = idx_start; j < idx_end; j++) {
                                depths[j]++;
                            }
                        }
                        
                        pos += op_len;
                        if (truncate_pos > 0 && pos >= truncate_pos) break;
                    } else if (cigar_consumes_ref[op]) {
                        pos += op_len;
                    }
                }
            }
        }
    }
    
    hts_itr_destroy(iter);
    
    // Compute stats directly
    uint64_t covered = 0;
    int64_t total = 0;
    const int32_t min_depth = args->min_depth;
    
    for (uint64_t i = 0; i < win_len; i++) {
        int32_t d = depths[i];
        total += d;
        if (d >= min_depth) covered++;
    }
    
    *covered_out = covered;
    *total_out = total;
    
    if (args->gc_content && gc_data && gc_data[tid]) {
        calculate_gc_content(gc_data[tid], win_start, win_end, gc_out, at_out);
    }
    
    return 0;
}

int process_chromosome_batch(int chr_tid, worker_context_t *ctx, region_t *chr_regions, int region_count) {
    if (region_count == 0) return 0;
    
    // Build region array for sam_itr_regarray
    char *char_map = (char*)malloc(region_count * 128);
    char **region_array = (char**)malloc(region_count * sizeof(char*));
    if (!char_map || !region_array) {
        free(char_map);
        free(region_array);
        return -1;
    }
    
    char *char_ptr = char_map;
    const char *chr_name = ctx->hdr->target_name[chr_tid];
    uint64_t chr_len = ctx->hdr->target_len[chr_tid];
    
    for (int i = 0; i < region_count; i++) {
        uint64_t beg = chr_regions[i].start;
        uint64_t end = chr_regions[i].end;
        if (beg < 1) beg = 1;
        if (end > chr_len) end = chr_len;
        
        region_array[i] = char_ptr;
        char_ptr += sprintf(char_ptr, "%s:%lu-%lu", chr_name, beg, end) + 1;
    }
    
    // Create multi-region iterator
    hts_itr_t *iter = sam_itr_regarray(ctx->idx, ctx->hdr, region_array, region_count);
    if (!iter) {
        free(char_map);
        free(region_array);
        return -1;
    }
    
    // Allocate depth array for entire chromosome span
    uint64_t min_pos = chr_regions[0].start;
    uint64_t max_pos = chr_regions[region_count - 1].end;
    uint64_t span = max_pos - min_pos;
    
    int32_t *depths = (int32_t*)calloc(span, sizeof(int32_t));
    if (!depths) {
        hts_itr_destroy(iter);
        free(char_map);
        free(region_array);
        return -1;
    }
    
    const uint32_t flags = ctx->args->exclude_flags;
    const uint8_t min_mapq = ctx->args->min_mapq;
    bam1_t **rec_batch = ctx->rec_batch;
    
    // Single iteration over all regions
    while (true) {
        int batch_count = 0;
        for (int b = 0; b < 128; b++) {
            if (sam_itr_next(ctx->fp, iter, rec_batch[b]) < 0) break;
            batch_count++;
        }
        if (batch_count == 0) break;
        
        for (int b = 0; b < batch_count; b++) {
            bam1_t *rec = rec_batch[b];
            
            if ((rec->core.flag & flags) || (rec->core.flag & BAM_FUNMAP) || rec->core.qual < min_mapq) {
                continue;
            }
            
            const uint32_t *cigar = bam_get_cigar(rec);
            int32_t pos = rec->core.pos - min_pos;
            
            for (uint32_t i = 0; i < rec->core.n_cigar; i++) {
                const int op = bam_cigar_op(cigar[i]);
                const int op_len = bam_cigar_oplen(cigar[i]);
                
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    int32_t end = pos + op_len;
                    for (int32_t j = pos; j < end; j++) {
                        if (j >= 0 && j < (int32_t)span) depths[j]++;
                    }
                    pos = end;
                } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                    pos += op_len;
                }
            }
        }
    }
    
    hts_itr_destroy(iter);
    free(char_map);
    free(region_array);
    
    // Extract statistics for each region
    const int32_t min_depth = ctx->args->min_depth;
    for (int r = 0; r < region_count; r++) {
        uint64_t win_start = chr_regions[r].start - min_pos;
        uint64_t win_len = chr_regions[r].end - chr_regions[r].start;
        int task_idx = chr_regions[r].region_idx;
        
        uint64_t covered = 0;
        int64_t total = 0;
        
        for (uint64_t i = 0; i < win_len; i++) {
            int32_t d = depths[win_start + i];
            total += d;
            if (d >= min_depth) covered++;
        }
        
        ctx->local_covered_sites[task_idx] += covered;
        ctx->local_total_depth[task_idx] += total;
        
        if (ctx->args->gc_content && ctx->gc_data && ctx->gc_data[chr_tid]) {
            uint64_t gc_bases = 0;
            uint64_t at_bases = 0;
            calculate_gc_content(ctx->gc_data[chr_tid], chr_regions[r].start, chr_regions[r].end, &gc_bases, &at_bases);
            ctx->local_gc_bases[task_idx] += gc_bases;
            ctx->local_at_bases[task_idx] += at_bases;
        }
    }
    
    return 0;
}

void* worker_thread_func(void *arg) {
    worker_context_t *ctx = (worker_context_t*)arg;
    region_t task;
    
    // Simple atomic scheduling
    while (atomic_scheduler_get_task(ctx->scheduler, &task)) {
        process_region(&task, ctx);
    }
    
    // Flush thread-local stats to global
    for (size_t i = 0; i < ctx->n_regions; i++) {
        if (ctx->local_covered_sites[i] > 0) {
            atomic_fetch_add_explicit(&ctx->region_stats[i].covered_sites,
                                     ctx->local_covered_sites[i], memory_order_relaxed);
        }
        if (ctx->local_total_depth[i] > 0) {
            atomic_fetch_add_explicit(&ctx->region_stats[i].total_depth,
                                     ctx->local_total_depth[i], memory_order_relaxed);
        }
        if (ctx->local_gc_bases[i] > 0) {
            atomic_fetch_add_explicit(&ctx->region_stats[i].gc_bases,
                                     ctx->local_gc_bases[i], memory_order_relaxed);
        }
        if (ctx->local_at_bases[i] > 0) {
            atomic_fetch_add_explicit(&ctx->region_stats[i].at_bases,
                                     ctx->local_at_bases[i], memory_order_relaxed);
        }
    }
    
    return NULL;
}

// ============================================================================
// MAIN
// ============================================================================

void print_banner() {
    printf(
        "\n"
        "   XX;          xx                                                                     \n"
        "   xXXXXXXXXXxxxxx                                                                     \n"
        "    XX;        +x;                                                                     \n"
        "     xXXXXxxxxxx;                                                                      \n"
        "      ;xXx;+xx+                                                                        \n"
        "        +xxxx           xX;                                                            \n"
        "      ;xx+:;xx+         X$+                                                            \n"
        "     xxxxxxxxxxx;       X$x$$$x ;$$    X$ X$xX$$X;   ;$$$$;  +$xX$; ;$$$$+   +$$$X  X$+   x$x  \n"
        "    xx;        +x;      X$X;x$$ x$X  x$x X$$; ;$$+ X$x  +$X +$$+; X$x  +X  $$+; X$x $$  ;$X   \n"
        "   xxxxxxxxx+++xx+      X$+   $$  X$x+$X  X$x   +$$ $$$$$$$$ +$x  ;$$      +$x   ;$$ ;$$ $$    \n"
        "   xx;                  X$+   $$  ;$$$$;  X$$+ ;$$+ x$X  +X  +$x   X$x  x$;;$$+; X$x  x$$$+    \n"
        "   ;xx      ;;;         xX;   XX   ;$$+   X$xx$$x     x$$$;  +Xx    ;X$$X;   ;$$$x     XXx     \n"
        "    ;xx     xxx                   +$$+    X$x                                                  \n"
        "      xx+   xxx +xx               xx;     +x;                                                  \n"
        "            xxx +xx                                                                            \n"
        "        xx+ xxx +xx                      hypercov - Hyper Fast Coverage Calculation            \n"
        "    ;;  xx+ xxx +xx                                                                            \n"
        "    xx; xx+ xxx +xx                                                                            \n"
        "\n"
    );
}

void print_usage(const char *prog) {
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  USAGE\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    fprintf(stderr, " %s -i <input.bam> -o <output_prefix> [options]\n\n", prog);
    
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  REQUIRED ARGUMENTS\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    fprintf(stderr, " -i <file>      Input alignment file (BAM/CRAM format)\n");
    fprintf(stderr, " -o <prefix>    Output file prefix for results\n\n");
    
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  ANALYSIS MODES (choose one, default = whole genome)\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    fprintf(stderr, " -b <file>      BED file - analyze specific genomic regions\n");
    fprintf(stderr, " -g <file>      GFF/GTF file - analyze gene features\n");
    fprintf(stderr, " -w <size>      Fixed windows - divide genome into windows of <size> bp\n");
    fprintf(stderr, " (default)      Whole genome - analyze entire genome coverage\n\n");
    
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  GENERAL OPTIONS\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    fprintf(stderr, " -r <file>      Reference FASTA (required for CRAM files)\n");
    fprintf(stderr, " -t <int>       Number of threads [default: 4]\n");
    fprintf(stderr, " -q <int>       Minimum mapping quality [default: 0]\n");
    fprintf(stderr, " -F <int>       Exclude reads with these flags [default: 1796]\n");
    fprintf(stderr, "                 (1796 = unmapped, secondary, QC fail, duplicate)\n");
    fprintf(stderr, " --min-depth <int>  Minimum depth to count as covered [default: 1]\n\n");
    
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  ADVANCED OPTIONS\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    fprintf(stderr, " -a             Generate per-base depth output (large file!)\n");
    fprintf(stderr, " -c             Calculate GC content (requires -r reference)\n");
    fprintf(stderr, " --unsorted     Process unsorted BAM (will sort on-the-fly)\n");
    fprintf(stderr, " --feature <type>   Feature type for GFF mode [default: gene]\n");
    fprintf(stderr, "                    Examples: gene, CDS, exon, mRNA\n\n");
    
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  OUTPUT FILES\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    fprintf(stderr, " <prefix>.coverage.gz       Main coverage statistics (gzip compressed)\n");
    fprintf(stderr, " <prefix>.per-base.bed.gz   Per-base depths (only with -a flag)\n\n");
    
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  EXAMPLES\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    fprintf(stderr, " # Whole genome coverage with 8 threads\n");
    fprintf(stderr, " %s -i sample.bam -o results -t 8\n\n", prog);
    
    fprintf(stderr, " # Analyze specific genomic regions from BED file\n");
    fprintf(stderr, " %s -i sample.bam -o results -b regions.bed -t 8\n\n", prog);
    
    fprintf(stderr, " # Gene coverage from GFF annotation\n");
    fprintf(stderr, " %s -i sample.bam -o results -g genes.gff -t 8\n\n", prog);
    
    fprintf(stderr, " # Coverage in 1Mb windows across genome\n");
    fprintf(stderr, " %s -i sample.bam -o results -w 1000000 -t 8\n\n", prog);
    
    fprintf(stderr, " # With GC content calculation\n");
    fprintf(stderr, " %s -i sample.bam -o results -r genome.fa -c -t 8\n\n", prog);
    
    fprintf(stderr, " # Process CRAM file with reference\n");
    fprintf(stderr, " %s -i sample.cram -r genome.fa -o results -t 8\n\n", prog);
    
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
}

int main(int argc, char **argv) {
    // Show banner
    print_banner();
    
    args_t args = {
        .input = NULL,
        .output = NULL,
        .reference = NULL,
        .bed_file = NULL,
        .gff_file = NULL,
        .threads = 4,
        .min_mapq = 0,
        .exclude_flags = 1796,
        .min_depth = 1,
        .window_size = 0,
        .all_sites = false,
        .unsorted = false,
        .gc_content = false,
        .mode = MODE_WHOLE_GENOME
    };
    
    char *feature_type = "gene";
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) args.input = argv[++i];
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) args.output = argv[++i];
        else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) args.reference = argv[++i];
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            args.bed_file = argv[++i];
            args.mode = MODE_BED_REGIONS;
        }
        else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            args.gff_file = argv[++i];
            args.mode = MODE_GFF_GENES;
        }
        else if (strcmp(argv[i], "-w") == 0 && i + 1 < argc) {
            args.window_size = atoll(argv[++i]);
            args.mode = MODE_FIXED_WINDOWS;
        }
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) args.threads = atoi(argv[++i]);
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) args.min_mapq = atoi(argv[++i]);
        else if (strcmp(argv[i], "-F") == 0 && i + 1 < argc) args.exclude_flags = atoi(argv[++i]);
        else if (strcmp(argv[i], "--min-depth") == 0 && i + 1 < argc) args.min_depth = atoi(argv[++i]);
        else if (strcmp(argv[i], "--feature") == 0 && i + 1 < argc) feature_type = argv[++i];
        else if (strcmp(argv[i], "-a") == 0) args.all_sites = true;
        else if (strcmp(argv[i], "-c") == 0) args.gc_content = true;
        else if (strcmp(argv[i], "--unsorted") == 0) args.unsorted = true;
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    if (!args.input || !args.output) {
        fprintf(stderr, "ERROR: Missing required arguments\n\n");
        print_usage(argv[0]);
        return 1;
    }
    
    double total_start = get_time();
    
    char *input_file = args.input;
    char *temp_sorted = NULL;
    
    if (args.unsorted) {
        temp_sorted = sort_bam_file(args.input, args.reference, args.threads);
        if (!temp_sorted) {
            fprintf(stderr, "Error: Failed to sort BAM file\n");
            return 1;
        }
        input_file = temp_sorted;
    }
    
    samFile *fp = sam_open(input_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", input_file);
        if (temp_sorted) {
            unlink(temp_sorted);
            free(temp_sorted);
        }
        return 1;
    }
    
    bool is_cram = is_cram_file(input_file);
    if (is_cram) {
        configure_input_file(fp, true, args.reference);
    }
    
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        fprintf(stderr, "Error: Cannot read header\n");
        sam_close(fp);
        if (temp_sorted) {
            unlink(temp_sorted);
            free(temp_sorted);
        }
        return 1;
    }
    
    sam_close(fp);
    
    region_list_t *regions = NULL;
    
    switch (args.mode) {
        case MODE_BED_REGIONS:
            regions = parse_bed_file(args.bed_file, hdr);
            break;
        case MODE_GFF_GENES:
            regions = parse_gff_file(args.gff_file, hdr, feature_type);
            break;
        case MODE_FIXED_WINDOWS:
            regions = generate_fixed_windows(hdr, args.window_size);
            break;
        case MODE_WHOLE_GENOME:
        default:
            regions = generate_processing_windows(hdr);
            break;
    }
    
    if (!regions || regions->count == 0) {
        fprintf(stderr, "Error: No regions to process\n");
        bam_hdr_destroy(hdr);
        if (temp_sorted) {
            unlink(temp_sorted);
            free(temp_sorted);
        }
        return 1;
    }
    
    uint64_t max_window_len = 0;
    for (size_t i = 0; i < regions->count; i++) {
        uint64_t len = regions->regions[i].end - regions->regions[i].start;
        if (len > max_window_len) max_window_len = len;
    }
    if (max_window_len == 0) max_window_len = 1;
    
    size_t stats_count = (args.mode == MODE_WHOLE_GENOME) ? (size_t)hdr->n_targets : regions->count;
    region_stats_t *region_stats = (region_stats_t*)calloc(stats_count, sizeof(region_stats_t));
    
    if (args.mode == MODE_WHOLE_GENOME) {
        for (int i = 0; i < hdr->n_targets; i++) {
            region_stats[i].name = strdup(hdr->target_name[i]);
            region_stats[i].total_bases = hdr->target_len[i];
            atomic_init(&region_stats[i].covered_sites, 0);
            atomic_init(&region_stats[i].total_depth, 0);
            atomic_init(&region_stats[i].gc_bases, 0);
            atomic_init(&region_stats[i].at_bases, 0);
        }
    } else {
        for (size_t i = 0; i < regions->count; i++) {
            region_stats[i].name = strdup(regions->regions[i].name ? regions->regions[i].name : "unknown");
            region_stats[i].total_bases = regions->regions[i].end - regions->regions[i].start;
            atomic_init(&region_stats[i].covered_sites, 0);
            atomic_init(&region_stats[i].total_depth, 0);
            atomic_init(&region_stats[i].gc_bases, 0);
            atomic_init(&region_stats[i].at_bases, 0);
        }
    }
    
    uint8_t **gc_data = NULL;
    if (args.gc_content) {
        if (!args.reference) {
            fprintf(stderr, "Error: GC content calculation requires -r reference FASTA\n");
            free(region_stats);
            free_region_list(regions);
            bam_hdr_destroy(hdr);
            if (temp_sorted) {
                unlink(temp_sorted);
                free(temp_sorted);
            }
            return 1;
        }
        gc_data = load_gc_content(args.reference, hdr);
        if (!gc_data) {
            fprintf(stderr, "Warning: Failed to load GC content, continuing without it\n");
            args.gc_content = false;
        }
    }
    
    char output_path[1024];
    snprintf(output_path, sizeof(output_path), "%s.coverage.gz", args.output);
    writer_thread_t *stats_writer = writer_thread_create(output_path);
    
    writer_msg_t header_msg = {
        .type = MSG_BATCH_LINES,
        .lines = (char**)malloc(sizeof(char*)),
        .count = 1
    };
    
    if (args.mode == MODE_GFF_GENES || args.mode == MODE_BED_REGIONS) {
        if (args.gc_content) {
            header_msg.lines[0] = strdup("#Region\tChrom\tStart\tEnd\tLength\tCoveredBases\tTotalDepth\tCoverage(%)\tMeanDepth\tGC(%)");
        } else {
            header_msg.lines[0] = strdup("#Region\tChrom\tStart\tEnd\tLength\tCoveredBases\tTotalDepth\tCoverage(%)\tMeanDepth");
        }
    } else {
        if (args.gc_content) {
            header_msg.lines[0] = strdup("#Chrom\tLength\tCoveredBases\tTotalDepth\tCoverage(%)\tMeanDepth\tGC(%)");
        } else {
            header_msg.lines[0] = strdup("#Chrom\tLength\tCoveredBases\tTotalDepth\tCoverage(%)\tMeanDepth");
        }
    }
    writer_thread_send(stats_writer, header_msg);
    
    writer_thread_t *perbase_writer = NULL;
    if (args.all_sites) {
        snprintf(output_path, sizeof(output_path), "%s.per-base.bed.gz", args.output);
        perbase_writer = writer_thread_create(output_path);
    }
    
    // OPTIMIZE: Single-threaded fast path
    if (args.threads == 1) {
        samFile *fp = sam_open(input_file, "r");
        if (!fp) { fprintf(stderr, "Error: Cannot open %s\n", input_file); return 1; }
        hts_set_threads(fp, 4);
        if (is_cram) {
            hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
            hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0);
            hts_set_opt(fp, HTS_OPT_CACHE_SIZE, CRAM_CACHE_SIZE);
            hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, 256 * 1024);
            if (args.reference) hts_set_fai_filename(fp, args.reference);
        }
        bam_hdr_t *hdr_local = sam_hdr_read(fp);
        hts_idx_t *idx_local = sam_index_load(fp, input_file);
        bam1_t **rec_batch = (bam1_t**)malloc(SINGLE_THREAD_BATCH_SIZE * sizeof(bam1_t*));
        for (int j = 0; j < SINGLE_THREAD_BATCH_SIZE; j++) rec_batch[j] = bam_init1();
        int32_t *depths = (int32_t*)aligned_alloc(CACHE_LINE_SIZE, max_window_len * sizeof(int32_t));
        
        // Local accumulators (non-atomic)
        uint64_t *local_cov = (uint64_t*)calloc(stats_count, sizeof(uint64_t));
        int64_t  *local_tot = (int64_t*)calloc(stats_count, sizeof(int64_t));
        uint64_t *local_gc  = (uint64_t*)calloc(stats_count, sizeof(uint64_t));
        uint64_t *local_at  = (uint64_t*)calloc(stats_count, sizeof(uint64_t));
        
        for (size_t r = 0; r < regions->count; r++) {
            region_t *task = &regions->regions[r];
            uint64_t covered = 0; int64_t total = 0; uint64_t gc = 0, at = 0;
            
            int ret = process_region_single_thread(task, fp, hdr_local, idx_local, rec_batch,
                                         depths, &args, &covered, &total, &gc, &at, gc_data);
            
            if (ret < 0) {
                fprintf(stderr, "Warning: Failed to process region %zu\n", r);
                continue;
            }
            
            // For whole genome mode, use tid to aggregate by chromosome
            // For other modes, use the linear index r
            int stats_idx = (args.mode == MODE_WHOLE_GENOME) ? task->tid : r;
            
            local_cov[stats_idx] += covered;
            local_tot[stats_idx] += total;
            if (args.gc_content) { local_gc[stats_idx] += gc; local_at[stats_idx] += at; }
            if ((r + 1) % 100 == 0) fprintf(stderr, "\rProcessed %zu/%zu regions", r + 1, regions->count);
        }
        fprintf(stderr, "\rProcessed %zu/%zu regions\n", regions->count, regions->count);
        
        // Final write-back
        for (size_t i = 0; i < stats_count; i++) {
            atomic_store_explicit(&region_stats[i].covered_sites, local_cov[i], memory_order_release);
            atomic_store_explicit(&region_stats[i].total_depth, local_tot[i], memory_order_release);
            if (args.gc_content) {
                atomic_store_explicit(&region_stats[i].gc_bases, local_gc[i], memory_order_release);
                atomic_store_explicit(&region_stats[i].at_bases, local_at[i], memory_order_release);
            }
        }
        
        free(local_cov); free(local_tot); free(local_gc); free(local_at);
        free(depths);
        for (int j = 0; j < SINGLE_THREAD_BATCH_SIZE; j++) bam_destroy1(rec_batch[j]);
        free(rec_batch);
        hts_idx_destroy(idx_local);
        bam_hdr_destroy(hdr_local);
        sam_close(fp);
    } else {
        // MULTI-THREADED PATH - keep existing code
        
        // Create atomic scheduler instead of work queue
        atomic_scheduler_t *scheduler = atomic_scheduler_create(regions->regions, regions->count);
        if (!scheduler) {
            fprintf(stderr, "Error: Cannot create scheduler\n");
            return 1;
        }
        
        worker_context_t *contexts = (worker_context_t*)calloc(args.threads, sizeof(worker_context_t));
        pthread_t *threads = (pthread_t*)malloc(args.threads * sizeof(pthread_t));
        
        for (int i = 0; i < args.threads; i++) {
            contexts[i].thread_id = i;
            contexts[i].fp = sam_open(input_file, "r");
            if (!contexts[i].fp) {
                fprintf(stderr, "Error: Thread %d cannot open input file\n", i);
                return 1;
            }
            
            contexts[i].is_cram = is_cram;
            
            hts_set_threads(contexts[i].fp, 3);
            
            if (is_cram) {
                hts_set_opt(contexts[i].fp, CRAM_OPT_REQUIRED_FIELDS,
                           SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
                hts_set_opt(contexts[i].fp, CRAM_OPT_DECODE_MD, 0);
                hts_set_opt(contexts[i].fp, HTS_OPT_CACHE_SIZE, CRAM_CACHE_SIZE);
                hts_set_opt(contexts[i].fp, HTS_OPT_BLOCK_SIZE, 256 * 1024);
                if (args.reference) {
                    hts_set_fai_filename(contexts[i].fp, args.reference);
                }
            } else {
                hts_set_opt(contexts[i].fp, HTS_OPT_BLOCK_SIZE, 256 * 1024);
                if (args.reference) {
                    hts_set_fai_filename(contexts[i].fp, args.reference);
                }
            }
            
            contexts[i].hdr = sam_hdr_read(contexts[i].fp);
            contexts[i].idx = sam_index_load(contexts[i].fp, input_file);
            
            contexts[i].rec_batch = (bam1_t**)malloc(READ_BATCH_SIZE * sizeof(bam1_t*));
            for (int j = 0; j < READ_BATCH_SIZE; j++) {
                contexts[i].rec_batch[j] = bam_init1();
            }
            
            contexts[i].buffer = simple_buffer_create(max_window_len);
            if (!contexts[i].buffer) {
                fprintf(stderr, "Error: Cannot allocate buffer for thread %d\n", i);
                return 1;
            }
            
            contexts[i].scheduler = scheduler;
            contexts[i].args = &args;
            contexts[i].region_stats = region_stats;
            contexts[i].n_regions = stats_count;
            contexts[i].perbase_writer = perbase_writer;
            contexts[i].gc_data = gc_data;
            
            contexts[i].local_covered_sites = (uint64_t*)calloc(stats_count, sizeof(uint64_t));
            contexts[i].local_total_depth = (int64_t*)calloc(stats_count, sizeof(int64_t));
            contexts[i].local_gc_bases = (uint64_t*)calloc(stats_count, sizeof(uint64_t));
            contexts[i].local_at_bases = (uint64_t*)calloc(stats_count, sizeof(uint64_t));
            
            contexts[i].local_reads_processed = 0;
            contexts[i].local_reads_skipped = 0;
        }
        
        for (int i = 0; i < args.threads; i++) {
            pthread_create(&threads[i], NULL, worker_thread_func, &contexts[i]);
        }
        
        uint64_t total_reads_processed = 0;
        uint64_t total_reads_skipped = 0;
        
        for (int i = 0; i < args.threads; i++) {
            pthread_join(threads[i], NULL);
            total_reads_processed += contexts[i].local_reads_processed;
            total_reads_skipped += contexts[i].local_reads_skipped;
        }
        
        // Cleanup multi-thread resources
        for (int i = 0; i < args.threads; i++) {
            for (int j = 0; j < READ_BATCH_SIZE; j++) {
                bam_destroy1(contexts[i].rec_batch[j]);
            }
            free(contexts[i].rec_batch);
            simple_buffer_free(contexts[i].buffer);
            free(contexts[i].local_covered_sites);
            free(contexts[i].local_total_depth);
            free(contexts[i].local_gc_bases);
            free(contexts[i].local_at_bases);
            hts_idx_destroy(contexts[i].idx);
            bam_hdr_destroy(contexts[i].hdr);
            sam_close(contexts[i].fp);
        }
        
        free(contexts);
        free(threads);
        atomic_scheduler_free(scheduler);
    }
    
    double total_elapsed = get_time() - total_start;
    
    size_t output_count = (args.mode == MODE_WHOLE_GENOME) ? (size_t)hdr->n_targets : regions->count;
    
    for (size_t i = 0; i < output_count; i++) {
        uint64_t total_bases = region_stats[i].total_bases;
        uint64_t covered = atomic_load_explicit(&region_stats[i].covered_sites, memory_order_acquire);
        int64_t total = atomic_load_explicit(&region_stats[i].total_depth, memory_order_acquire);
        uint64_t gc = atomic_load_explicit(&region_stats[i].gc_bases, memory_order_acquire);
        uint64_t at = atomic_load_explicit(&region_stats[i].at_bases, memory_order_acquire);
        
        double pct = (total_bases > 0) ? ((double)covered / (double)total_bases) * 100.0 : 0.0;
        double mean = (total_bases > 0) ? ((double)total / (double)total_bases) : 0.0;
        double gc_pct = (gc + at > 0) ? ((double)gc / (double)(gc + at)) * 100.0 : 0.0;
        
        char line[1024];
        if (args.mode == MODE_WHOLE_GENOME) {
            if (args.gc_content) {
                snprintf(line, sizeof(line), "%s\t%lu\t%lu\t%ld\t%.2f\t%.2f\t%.2f",
                        region_stats[i].name, total_bases, covered, total, pct, mean, gc_pct);
            } else {
                snprintf(line, sizeof(line), "%s\t%lu\t%lu\t%ld\t%.2f\t%.2f",
                        region_stats[i].name, total_bases, covered, total, pct, mean);
            }
        } else if (args.mode == MODE_GFF_GENES || args.mode == MODE_BED_REGIONS) {
            if (args.gc_content) {
                snprintf(line, sizeof(line), "%s\t%s\t%lu\t%lu\t%lu\t%lu\t%ld\t%.2f\t%.2f\t%.2f",
                        region_stats[i].name, regions->regions[i].chrom_name,
                        regions->regions[i].start, regions->regions[i].end,
                        total_bases, covered, total, pct, mean, gc_pct);
            } else {
                snprintf(line, sizeof(line), "%s\t%s\t%lu\t%lu\t%lu\t%lu\t%ld\t%.2f\t%.2f",
                        region_stats[i].name, regions->regions[i].chrom_name,
                        regions->regions[i].start, regions->regions[i].end,
                        total_bases, covered, total, pct, mean);
            }
        } else {
            if (args.gc_content) {
                snprintf(line, sizeof(line), "%s\t%lu\t%lu\t%lu\t%lu\t%ld\t%.2f\t%.2f\t%.2f",
                        regions->regions[i].chrom_name,
                        regions->regions[i].start, regions->regions[i].end,
                        total_bases, covered, total, pct, mean, gc_pct);
            } else {
                snprintf(line, sizeof(line), "%s\t%lu\t%lu\t%lu\t%lu\t%ld\t%.2f\t%.2f",
                        regions->regions[i].chrom_name,
                        regions->regions[i].start, regions->regions[i].end,
                        total_bases, covered, total, pct, mean);
            }
        }
        
        writer_msg_t msg = {
            .type = MSG_BATCH_LINES,
            .lines = (char**)malloc(sizeof(char*)),
            .count = 1
        };
        msg.lines[0] = strdup(line);
        writer_thread_send(stats_writer, msg);
    }
    
    writer_thread_close(stats_writer);
    if (perbase_writer) writer_thread_close(perbase_writer);
    
    // Final cleanup
    for (size_t i = 0; i < stats_count; i++) {
        free(region_stats[i].name);
    }
    
    free(region_stats);
    free_region_list(regions);
    if (gc_data) free_gc_data(gc_data, hdr->n_targets);
    bam_hdr_destroy(hdr);
    
    if (temp_sorted) {
        char idx_file[512];
        snprintf(idx_file, sizeof(idx_file), "%s.bai", temp_sorted);
        unlink(idx_file);
        unlink(temp_sorted);
        free(temp_sorted);
    }
    
    fprintf(stderr, "\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  ✓ Coverage calculation completed successfully!\n");
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n");
    fprintf(stderr, "  Time elapsed: %.2f seconds\n", total_elapsed);
    fprintf(stderr, "  Output files:\n");
    fprintf(stderr, "    - %s.coverage.gz\n", args.output);
    if (args.all_sites) {
        fprintf(stderr, "    - %s.per-base.bed.gz\n", args.output);
    }
    fprintf(stderr, "═══════════════════════════════════════════════════════════════════════\n\n");
    
    return 0;
}