all: fqzcomp5

# Compiler and flags
CC ?= gcc
INCLUDES = -I.
CFLAGS = -g -O3 -Wall
LDFLAGS =
LIBS = -lz -lm -lbz2 -pthread

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

clean:
	-rm -f fqzcomp5 $(FQZCOMP5_OBJ) $(ALL_HTSCODECS_OBJ)

.PHONY: all clean
