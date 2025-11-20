<p align="left">
  <img src="./logo hypercov.png" alt="HyperCov Logo" width="75%" />
</p>

HyperCov is a high-performance coverage calculator. It supports sorted and unsorted BAM files as well as CRAM files. Hypercov provides several coverage calculation modes supporting BED/GFF regions, fixed windows, whole-genome mode, GC calculation and per-base output.

## Installation

Clone the repository and build from source, or use a provided precompiled binary.

- Clone:
```bash
git clone https://github.com/lisosome/hypercov.git
cd hypercov
```

- Expected local dependency layout (Makefile relies on these paths):
  - Headers: `libs/include/` (e.g. `libs/include/sam.h`, `libs/include/hts.h`, `libs/include/libdeflate.h`)
  - Libraries: `libs/lib/` (e.g. `libs/lib/libhts.a` or `libhts.so`, `libs/lib/libdeflate.a` or `libdeflate.so`)

- Build:
```bash
make
```

- Precompiled binary:
  - A precompiled `hypercov` binary may be provided in repository Releases or placed in the repo root. If present you can run it directly:

```bash
./hypercov --help
```

## Usage
Basic:
```bash
./hypercov -i sample.bam -o results
```

Common options
- `-i <file>`    Input alignment file (BAM/CRAM)
- `-o <prefix>`  Output file prefix (creates `<prefix>.coverage.gz`)
- `-r <file>`    Reference FASTA (required for CRAM and GC)
- `-t <int>`     Threads (default: 4)
- `-q <int>`     Minimum mapping quality
- `-F <int>`     Exclude reads with these flags (default: 1796)
- `-b <file>`    BED file (region list)
- `-g <file>`    GFF/GTF file (gene/features)
- `-w <size>`    Fixed window size (bp)
- `-a`           Produce per-base output (large)
- `-c`           Compute GC content (requires `-r`)
- `--unsorted`   Sort BAM on-the-fly if not sorted
- `--min-depth <int>` minimum depth to consider a base covered (default: 1)

## Examples
Whole-genome with 8 threads:
```bash
./hypercov -i sample.bam -o results -t 8
```
BED regions:
```bash
./hypercov -i sample.bam -o results -b regions.bed -t 8
```
GFF gene coverage:
```bash
./hypercov -i sample.bam -o results -g genes.gff -t 8
```
1Mb windows:
```bash
./hypercov -i sample.bam -o results -w 1000000 -t 8
```
With GC content (requires reference):
```bash
./hypercov -i sample.bam -o results -r genome.fa -c -t 8
```

## Output
- `<prefix>.coverage.gz` — main coverage summary (gzip)
- `<prefix>.per-base.bed.gz` — per-base depths (when `-a`)
