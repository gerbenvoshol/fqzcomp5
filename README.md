Fqzcomp5 (ALPHA)
========

This is a reimplementation of the old fqzcomp-4.x tool using the new
htscodecs library (with a few amendments, to be folded back upstream).

Given the parameterisation of the fqzcomp qual model, this is
currently a little slower than the v4.6 code, but this may potentially
be fixable by compiling custom implementations for regularly used
configurations.

The file format is block based, typically between 100MB and
1GB in size.  This permits it to be more aggressively multi-threaded,
but it harms the maximum compression ratio slightly.

**NEW**: The file format now includes:
- A proper magic number (`FQZ5`) and version identifier
- An index at the end of the file for random access capability
- Better block structure with explicit block size metadata

This permits coarse random access, which for FASTQ files practically
means sub-dividing the data into blocks, e.g. for parallel dispatching
to an aligner. The index enables seeking directly to specific blocks
without reading the entire file.

Compared to fqzcomp-4:

- New fast modes (-1 and -3) using basic bit-packing and rANS.  This
  is designed for rapid fastq compression, and is also appropriate for
  smaller block sizes.

- More flexible fqzcomp_qual modelling, which on some data sets can
  give substantially smaller quality compression depending on the
  nature of the data.

- The optional addition of sequence based in the quality compression
  context.  This improves compression ratios for ONT and PacBio data
  sets.

- A learning algorithm (similar to our CRAM implementation) which
  performs periodic trials to evaluate what codecs work best.  The
  choice of codecs permitted is one of the things tweaked by -1 to -9
  (only odd numbers implemented at the moment).  For example -9 has
  the most number of alternatives for fqzcomp_qual models.

Building and Installation
==========================

Dependencies
------------

Fqzcomp5 requires the following dependencies:

- **zlib** - Compression library for gzipped FASTQ support
- **libbz2** - BZip2 compression library
- **gcc** or compatible C compiler

The htscodecs library is now fully integrated into the fqzcomp5 source tree,
so no separate installation or submodule initialization is required.

On Debian/Ubuntu systems, install dependencies with:
```bash
sudo apt-get install build-essential zlib1g-dev libbz2-dev
```

On Red Hat/Fedora systems:
```bash
sudo dnf install gcc zlib-devel bzip2-devel
```

Or for older systems using yum:
```bash
sudo yum install gcc zlib-devel bzip2-devel
```

Building from Source
--------------------

1. **Clone the repository:**
   ```bash
   git clone https://github.com/gerbenvoshol/fqzcomp5.git
   cd fqzcomp5
   ```

2. **Build the project:**
   ```bash
   make
   ```

   The Makefile will automatically:
   - Compile the integrated htscodecs library sources
   - Compile fqzcomp5 and link it with htscodecs
   - Create the `fqzcomp5` binary in the current directory

3. **Run tests:**
   ```bash
   ./test.sh
   ```

About htscodecs Integration
----------------------------

Fqzcomp5 uses the [htscodecs library](https://github.com/jkbonfield/htscodecs), which is now 
fully integrated into the source tree. This library provides the core compression codecs including:

- **rANS** (Asymmetric Numeral Systems) entropy encoders - both 4x8 and 4x16 variants with bit-packing/RLE
- **Adaptive arithmetic coding** for high-compression scenarios
- **FQZComp quality compression** - specialized quality score compression
- **Name tokenizer** - efficient read name compression
- **SIMD optimizations** - SSE4, AVX2, and AVX512 support for faster encoding/decoding

The integrated htscodecs code is based on the `fqz_seq_u32` branch from the upstream repository, 
which includes enhancements for sequence-based quality compression contexts. These enhancements 
significantly improve compression ratios for ONT and PacBio data.

The htscodecs sources are located in the `htscodecs/` directory and are compiled and statically 
linked into fqzcomp5 during the build process. No separate installation of htscodecs is required.

Usage
=====

Basic Usage
-----------

Fqzcomp5 supports compression and decompression of FASTQ files, with automatic handling of 
gzipped input/output and paired-end files.

**Compress a single FASTQ file:**
```bash
./fqzcomp5 input.fastq output.fqz5
```

**Decompress:**
```bash
./fqzcomp5 -d input.fqz5 output.fastq
```

**Compress paired-end files (with automatic interleaving):**
```bash
./fqzcomp5 input_R1.fastq input_R2.fastq output.fqz5
```

**Decompress to paired-end files:**
```bash
./fqzcomp5 -d input.fqz5 output_R1.fastq output_R2.fastq
```

**Gzipped files are handled transparently:**
```bash
# Compress gzipped FASTQ
./fqzcomp5 input.fastq.gz output.fqz5

# Decompress to gzipped FASTQ
./fqzcomp5 -d input.fqz5 output.fastq.gz

# Paired-end with gzipped files
./fqzcomp5 input_R1.fastq.gz input_R2.fastq.gz output.fqz5
```

**Using stdin/stdout:**
```bash
# Compress from stdin
cat input.fastq | ./fqzcomp5 > output.fqz5

# Decompress to stdout
./fqzcomp5 -d input.fqz5 | less

# Pipeline example
zcat input.fastq.gz | ./fqzcomp5 | ./fqzcomp5 -d > output.fastq
```

Compression Levels
------------------

Use `-1` through `-9` to control compression level (odd numbers only):

- **`-1`**: Fastest compression using basic bit-packing and rANS
- **`-3`**: Fast compression with rANS/LZP
- **`-5`**: Balanced compression with Markov models (default)
- **`-7`**: High compression with adaptive models
- **`-9`**: Maximum compression (slowest)

```bash
# Fast compression for quick archival
./fqzcomp5 -1 input.fastq output.fqz5

# Maximum compression for long-term storage
./fqzcomp5 -9 input.fastq output.fqz5
```

Advanced Options
----------------

```bash
# Use more threads (default is 4)
./fqzcomp5 -t 8 input.fastq output.fqz5

# Custom block size (default varies by compression level)
./fqzcomp5 -b 500M input.fastq output.fqz5

# Include read name on third line (+name instead of +)
./fqzcomp5 -d -p input.fqz5 output.fastq

# Adjust verbosity
./fqzcomp5 -v input.fastq output.fqz5  # More verbose
./fqzcomp5 -V input.fastq output.fqz5  # Silent mode

# Fine-tune encoding methods (advanced users)
./fqzcomp5 -n 2 -s 1 -q 1 input.fastq output.fqz5
```

File Format
===========

The FQZ5 file format consists of:

1. **Header** (16 bytes):
   - Magic number: `FQZ5\001\000\000\000` (8 bytes) - identifies file as FQZ5 version 1.0.0
   - Index offset: 8-byte unsigned integer pointing to the index location (0 if no index)

2. **Data Blocks** (variable size):
   Each block contains:
   - Block size (4 bytes): Size of block data excluding this field
   - Number of records (4 bytes)
   - Compressed name data
   - Compressed length data
   - Compressed sequence data
   - Compressed quality data

3. **Index** (optional, at end of file):
   - Index magic: `FQZ5IDX\000` (8 bytes)
   - Number of blocks (4 bytes)
   - For each block:
     - File offset (8 bytes)
     - Uncompressed size in bases (4 bytes)
     - Number of records (4 bytes)

The index enables random access by allowing seekers to jump directly to
any block without decompressing previous blocks.

**Backward Compatibility**: Files without the FQZ5 magic header are
automatically detected and processed using the old format logic.

Results
=======

In the results below I mostly took the standard options for the tools.
It's likely better results are available by changing more parameters,
such as increasing the block sizes.  I did not explore this however
as it's already a large search space.


ERR174310
---------

Illumina HiSeq 2000 paired end sequencing
Each of the two files compressed independently and results summed together.

This file is also MPEG-G Dataset 01 and has figures reported in their
paper.

| Tool                 | Total(MB)  | Name | Lengths  | Sequence | Qual  |  Elapsed(s) | CPU (s) |
| -------------------- | ---------: | ---: | -------: | -------: | ----: | ----------: | ------: |
| gzip                 |      36652 |      |          |          |       |             |         |
| mpeg-g (no assembly) |      24112 |  376 | (in seq) |    10249 | 13486 |  (-@8) 1193 |    6166 |
| fqzcomp5 -1          |      26652 | 2743 |        0 |    10104 | 13813 |         904 |    1145 |
| fqzcomp5 -3          |      25682 | 1770 |        0 |    10103 | 13808 |        1111 |    1600 |
| fqzcomp5 -5          |      23355 | 1770 |        0 |     9548 | 12037 |        1542 |    5472 |
| fqzcomp5 -7          |      21922 |  668 |        0 |     9288 | 11966 |        2009 |    6576 |
| fqzcomp5 -9          |      21758 |  667 |        0 |     9136 | 11956 |        2435 |    8567 |
| fqzcomp4 n2 s7+b q3  |      21009 |  817 | (in seq) |     8432 | 11759 |        4901 |    6967 |
| genozip              |      23726 |~1426 |        0 |    ~9771 |~12455 |        5472 |   21620 |
| genozip --best=NO_REF|      22995 |~1426 |        0 |    ~9448 |~12241 |       10376 |   39483 |
| CRAM 3.0             |      25121 | 1100 |  (fixed) |    10139 | 13867 |         908 |    1988 |
| CRAM 3.1 small       |      23231 |  404 |  (fixed) |    10102 | 12719 |         964 |    4740 |
| Spring               |      15368 |  407 |  (fixed) |     2431 | 12523 |       10943 |   35438 |

Fqzcomp5's elapsed time is with 4 threads.  GenomSys publication uses 8
threads, with the CPU time being reported as their elapsed 1-thread
benchmark as an estimation.  Both reports used very similar CPUs, but
the I/O systems may be different.  It's clear fqzcomp was I/O bound
for faster methods, and possibly their tool too.

Fqzcomp4 is the smallest, but it lacks random access capability.

The faster rANS/LZP based fqzcomp5 levels (-1 and -3) are still
reasonably competitive and clearly beat the original two fastq.gz
sizes (reported in gzip row).  Once adaptive range coding models are
used the CPU time is considerably slower, but data size drops
considerably.  Level 9 doesn't seem to offer any benefit here over 7.

Genozip-13.0.5 was ran in the default mode (as "best" needs a
reference).  Sizes are approximate due to coarse rounding in the
genocat --stats output.

CRAM sizes are a little different to fqzcomp here as the command used was

```
   samtools import -N -1 ERR174310_1.fastq.gz -2 ERR174310_2.fastq.gz
```

This interleaves the two files rather than compressing each separately
and aggregating results, which in turn halves the size of the read
name encoding due to deduplication.  It also performs reasonably well
with CRAM 3.1, but has poorer quality compression.  This is perhaps
due to finer-grained random access (unverified).   CRAM 3.1 lacks an
adaptive sequence model, so is behind fqzcomp on sequence
compression.

Spring also interleaves both files.  Spring is doing a local read
clustering and reordering process to help compress the sequence data.
This is both CPU and memory intensive (about 21GB on this data set).
While it can achieve good results, I would argue that the resources
would be best spent on sequence alignment and/or a proper denovo
assembly (or maybe a mixture of both, with assembly on the data that
doesn't map to identify the large insertions, organelles and
contaminants).  This is probably a significant step up again in
resources, but the end result is far more useful to the user.

I view FASTQ largely as an interim data format - somewhere between the
raw instrument data and the end product (assembled BAM/CRAM and/or
VCFs).  It doesn't make a great deal of sense to me to spend a huge
amount of CPU on FASTQ compression alone.


SRR1238539
----------

IonTorrent variable length WGS data.

This file is also MPEG-G Dataset 11 and has figures reported in their
paper.


| Tool                 | Total(MB)  | Name | Lengths  | Sequence | Qual  |  Elapsed(s) | CPU (s) |
| -------------------- | ---------: | ---: | -------: | -------: | ----: | ----------: | ------: |
| gzip                 |      25000 |      |          |          |       |             |    5608 |
| zstd -1              |      26590 |      |          |          |       |             |    2338 |
| zstd -6              |      24110 |      |          |          |       |             |    3120 |
| zstd -15             |      23205 |      |          |          |       |             |   25038 |
| xz                   |      21324 |      |          |          |       |             |   57910 |
| mpeg-g (no assembly) |      21769 |  205 | (in seq) |     7873 | 13692 |         960 |    5982 |
| fqzcomp5 -1          |      19391 |  818 |      337 |     3742 | 14494 |         378 |     748 |
| fqzcomp5 -3          |      19145 |  583 |      337 |     3733 | 14492 |         441 |     941 |
| fqzcomp5 -5          |      15910 |  583 |      337 |     3733 | 11257 |         713 |    2978 |
| fqzcomp5 -7          |      15353 |  163 |      337 |     3737 | 11116 |         998 |    3484 |
| fqzcomp5 -9          |      15332 |  163 |      337 |     3738 | 11093 |        1063 |    3673 |
| fqzcomp4 n2 s7 b q3  |      20319 |  187 | (in seq) |     6728 | 13241 |        2075 |    3305 |
| colord -balanced     |      17087 |  401 | (in seq) |     4043 | 12643 |        5480 |   23461 |
| colord -ratio        |      16033 |  401 | (in seq) |     2712 | 12921 |       10642 |   37808 |
| genozip              |      18926 | ~480 | (in seq) |    ~5261 |~13207 |        3776 |   14682 |
| genozip --best=NO_REF|      18511 | ~460 | (in seq) |    ~5154 |~12885 |        6509 |   25288 |
| CRAM 3.0             |      18340 |  316 |      209 |     3320 | 14490 |         370 |    1719 |
| CRAM 3.1 small       |      17302 |   ?  |      183 |     3340 | 13774 |        1239 |    5169 |

As before Fqzcomp's elapsed time is with 4 threads and MPEG-G's is 8 threads.

Fqzcomp's sequence is small here due to the use of LZP and an apparent
genomic ordering in sequences despite being a FASTQ file.  This data
was originally converted by EBI from an SRA file.  Maybe it's possible
SRA has already done some sequence sorting for compression purposes,
but the original read identifiers have been lost by the SRA so it's
not possible to glean the original ordering to verify this.

The gzip size here is as reported by rerunning gzip rather than taking
the input fastq.gz size (it started as .sra anyway), which means we
have times.  Also shown are other generic compression tools.  Xz is
reasonably competitive with light-weight fqz encoding, but is
extremely slow.

Fqzcomp -1 and -3 are fast rANS / LZP / tokenisation strategies.  -5
onwards is using markov modelling and adaptive range coding for
quality values, which markedly reduces the size.  The quality encoding
here is exploiting both previous quality values and neighbouring
sequence bases as context.

Fqzcomp4 performs poorly on this data set, compressing neither the
sequence nor the quality well.

Colord in balanced mode does a reasonable job, but is let down mainly
by the quality compression and the slow speed.  It doesn't have an
IonTorrent profile, but testing on a smaller subset the ONT profile
gave the best ratio so this was used for the full dataset. (On the
first 1 million records, the ONT profile is 10% smaller than the
others for quality values.)  I'm sure it would be more competitive on
sequence compression had the data not apparently already been sorted.

Genozip-13.0.5 was ran in the default mode and "--best=NO_REF".

CRAM on this data set performs well at lighter levels in CRAM 3.0,
again due to the unusual nature of the data already being
clustered by sequence identity.  CRAM 3.1 currently lacks the ability
to use sequence bases in the FQZ-qual model, so loses out considerably
at higher compression levels.


ERR2442595
----------

This is a relatively small RNASeq dataset sequenced on Illumina
NovaSeq 6000.  The quality on this instrument is binned to just 4
discrete values.  The nature of the RNAseq experiment gives an
effective small portion of genome covered, which aids sequence
compression.

| Tool                 | Total(MB)  | Name | Lengths  | Sequence | Qual  |  Elapsed(s) | CPU (s) |
| -------------------- | ---------: | ---: | -------: | -------: | ----: | ----------: | ------: |
| gzip                 |       3852 |      |          |          |       |             |         |
| fqzcomp5 -1          |       2625 |  132 |        0 |     2076 |   417 |         118 |     151 |
| fqzcomp5 -3          |       2538 |   60 |        0 |     2067 |   411 |         121 |     204 |
| fqzcomp5 -5          |       1817 |   60 |        0 |     1361 |   396 |         236 |     873 |
| fqzcomp5 -7          |       1406 |    0 |        0 |     1014 |   392 |         407 |    1359 |
| fqzcomp5 -9          |       1193 |    0 |        0 |      801 |   392 |         696 |    2316 |
| fqzcomp4 n2 s7 b q3  |       1038 |    0 | (in seq) |      636 |   402 |         973 |    1022 |
| spring               |        645 |    0 |        0 |      233 |   410 |         682 |    1970 |
| genozip --best=NO_REF|       1464 |   60 |        0 |      988 |   416 |         702 |    1926 |


Spring is included here to show the impact of sequence reordering and
clustering.  It significantly reduces compressed sequence size over
the fqzcomp statistical model.  On this data set it's doesn't have a
large memory and CPU usage, although this can become an issue with
bigger datasets as seen above


TO DO
=====

- ~~Include support for input and output of gzipped FASTQ.~~

  **DONE**: gzipped FASTQ files are now supported transparently. The tool
  automatically detects `.gz` file extensions and handles compression/decompression
  accordingly using kseq.h and zlib.
  
  Examples:
  - Compress gzipped FASTQ: `fqzcomp5 input.fastq.gz output.fqz5`
  - Decompress to gzipped FASTQ: `fqzcomp5 -d input.fqz5 output.fastq.gz`
  - Works with stdin/stdout: `zcat input.fastq.gz | fqzcomp5 > output.fqz5`

- ~~Accept pairs of fastq files and do automatic interleaving.  This
  helps improve compression of read-names through deduplication and
  may also help sequence compression for short inserts.~~

  **DONE**: Paired-end FASTQ files are now supported with automatic interleaving
  during compression and deinterleaving during decompression. This significantly
  improves compression ratios for read names through deduplication.
  
  Examples:
  - Compress paired files: `fqzcomp5 input_R1.fastq input_R2.fastq output.fqz5`
  - Decompress to paired files: `fqzcomp5 -d input.fqz5 output_R1.fastq output_R2.fastq`
  - Also works with gzipped files: `fqzcomp5 input_R1.fastq.gz input_R2.fastq.gz output.fqz5`

- ~~Support the third line of "+name" where "name" is a duplicate of
  "@name".  This is rarely used, but could be supported by a simple
  flag.~~

  **DONE**: The `-p` flag enables outputting the read name on the third line
  (e.g., `+name` instead of just `+`). While this format variant is rarely used,
  it's now supported for compatibility.
  
  Example:
  - Decompress with names on third line: `fqzcomp5 -d -p input.fqz5 output.fastq`

- ~~Improve file format. A proper magic number, better blocking structure.~~

  **DONE**: The file format now includes:
  - Magic number `FQZ5` with version identifier (version 1.0.0)
  - Explicit block size metadata in each block
  - Index offset in the header for quick index lookup
  
  The new format maintains backward compatibility - files without the magic
  number are automatically detected and handled using the legacy format.

- ~~Implement (coarse) random access capability. This is already
  supported by the format, but lacks the necessary index.~~

  **DONE**: Random access is now implemented via an index structure:
  - Index is written at the end of each compressed file
  - Contains block offsets, sizes, and record counts
  - Enables seeking to specific blocks without decompressing preceding data
  - Useful for parallel processing or selective decompression
  
  The index structure includes:
  - File offset of each block
  - Uncompressed size (number of bases) per block
  - Number of records per block
  
  See the "File Format" section above for details.

- ~~Push more changes back into htscodecs upstream.  For now we're using
  our own fork.~~

  **DONE**: The htscodecs library is now fully integrated into the fqzcomp5 source tree.
  The integrated code is based on the `fqz_seq_u32` branch from the upstream repository 
  (https://github.com/jkbonfield/htscodecs), which includes:
  - Support for using sequence bases as quality compression context
  - SMALL_MODEL and SIMPLE_MODEL optimizations
  - All the necessary enhancements for fqzcomp5
  
  The htscodecs sources are located in the `htscodecs/` directory and are automatically
  compiled and statically linked during the build process. No submodule or separate 
  installation is required. See the "Building and Installation" section for details.

- Fuzz testing.
