all: fqzcomp5

# Compiler and flags
CC ?= gcc
INCLUDES = -I.
CFLAGS = -O3 -Wall -march=native -mtune=native -ffast-math -funroll-loops -fomit-frame-pointer
LDFLAGS = -Wl,-O1 -Wl,--as-needed
LIBS = -lz -lm -lbz2 -pthread

# Debug build flags
DEBUG_CFLAGS = -g -O0 -Wall
DEBUG_LDFLAGS =  # Empty: debug builds use no special linker flags

# Static build flags
STATIC_LDFLAGS = -static
STATIC_BUILD_LDFLAGS = -Wl,-O1 -Wl,--as-needed $(STATIC_LDFLAGS)

# Portable flags for static builds (run on a wide range of Intel/AMD x86-64 CPUs)
PORTABLE_CFLAGS = -march=x86-64 -mtune=generic

# Main fqzcomp5 objects
FQZCOMP5_OBJ = fqzcomp5.o lzp16e.o thread_pool.o

# htscodecs base objects
HTSCODECS_OBJ = \
	htscodecs/pack.o \
	htscodecs/rle.o \
	htscodecs/fqzcomp_qual.o \
	htscodecs/rANS_static.o \
	htscodecs/rANS_static4x16pr.o \
	htscodecs/rANS_static32x16pr.o \
	htscodecs/tokenise_name3.o \
	htscodecs/arith_dynamic.o \
	htscodecs/htscodecs.o \
	htscodecs/utils.o

# SIMD optimized objects (compiled with specific flags)
HTSCODECS_SIMD_OBJ = \
	htscodecs/rANS_static32x16pr_sse4.o \
	htscodecs/rANS_static32x16pr_avx2.o \
	htscodecs/rANS_static32x16pr_avx512.o

# All htscodecs objects
ALL_HTSCODECS_OBJ = $(HTSCODECS_OBJ) $(HTSCODECS_SIMD_OBJ)

# Build fqzcomp5
fqzcomp5: $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ)
	$(CC) $(LDFLAGS) $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ) -o $@ $(LIBS)

# Main fqzcomp5 source files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# htscodecs base objects
htscodecs/%.o: htscodecs/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# SIMD-specific compilation rules
htscodecs/rANS_static32x16pr_sse4.o: htscodecs/rANS_static32x16pr_sse4.c
	$(CC) $(CFLAGS) $(INCLUDES) -msse4.1 -mssse3 -mpopcnt -c $< -o $@

htscodecs/rANS_static32x16pr_avx2.o: htscodecs/rANS_static32x16pr_avx2.c
	$(CC) $(CFLAGS) $(INCLUDES) -mavx2 -c $< -o $@

htscodecs/rANS_static32x16pr_avx512.o: htscodecs/rANS_static32x16pr_avx512.c
	$(CC) $(CFLAGS) $(INCLUDES) -mavx512f -c $< -o $@

# Debug build (with debug symbols and no optimization)
# Note: Run 'make clean' first when switching between build types
debug: CFLAGS = $(DEBUG_CFLAGS)
debug: LDFLAGS = $(DEBUG_LDFLAGS)
debug: fqzcomp5-debug

fqzcomp5-debug: $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ)
	$(CC) $(LDFLAGS) $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ) -o $@ $(LIBS)

# Static build (statically linked executable)
# Note: Run 'make clean' first when switching between build types
static: CFLAGS = $(filter-out -march=native -mtune=native,$(CFLAGS)) $(PORTABLE_CFLAGS)
static: LDFLAGS = $(STATIC_BUILD_LDFLAGS)
static: fqzcomp5-static

fqzcomp5-static: $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ)
	$(CC) $(LDFLAGS) $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ) -o $@ $(LIBS)

# Debug static build (debug symbols + statically linked, no linker optimizations)
# Note: Run 'make clean' first when switching between build types
debug-static: CFLAGS = $(DEBUG_CFLAGS) $(PORTABLE_CFLAGS)
debug-static: LDFLAGS = $(DEBUG_LDFLAGS) $(STATIC_LDFLAGS)  # No -Wl,-O1 for easier debugging
debug-static: fqzcomp5-debug-static

fqzcomp5-debug-static: $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ)
	$(CC) $(LDFLAGS) $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ) -o $@ $(LIBS)

clean:
	-rm -f fqzcomp5 fqzcomp5-debug fqzcomp5-static fqzcomp5-debug-static $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ)

.PHONY: all clean debug static debug-static
