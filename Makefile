# Directories
SRC_DIR = src
LIBS_DIR = libs

# Compiler and flags
CC = gcc
CFLAGS = -O3 -march=native -mtune=native -mavx2 -mfma -flto=auto -ffast-math \
         -funroll-loops -finline-functions \
         -I$(LIBS_DIR)/include -I$(LIBS_DIR)/lib

LDFLAGS = -L$(LIBS_DIR)/lib \
          -Wl,-rpath,$(LIBS_DIR)/lib

LIBS = -lhts -ldeflate -lpthread -lm -lz

# Target
TARGET = hypercov
SRC = $(SRC_DIR)/hypercov.c

# Build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS) $(LIBS)
	@echo "Build complete: $(TARGET)"

# Install target (optional)
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/
	@echo "Installed to /usr/local/bin/$(TARGET)"

# Clean
clean:
	rm -f $(TARGET)

# Test
test: $(TARGET)
	@echo "Running basic test..."
	./$(TARGET) --help

.PHONY: all clean install test

