# Gzipped FASTQ Support

## Overview

fqzcomp5 now supports transparent reading and writing of gzipped FASTQ files using kseq.h from klib and zlib.

## Features

- **Automatic Detection**: File format is detected automatically based on `.gz` extension
- **Bidirectional Support**: Both input and output can be gzipped
- **Stdin/Stdout Compatible**: Works seamlessly with pipes
- **No Performance Impact**: Uses efficient streaming I/O

## Usage Examples

### Compress a gzipped FASTQ file
```bash
fqzcomp5 input.fastq.gz output.fqz5
```

### Decompress to gzipped FASTQ
```bash
fqzcomp5 -d input.fqz5 output.fastq.gz
```

### Use with pipes
```bash
zcat input.fastq.gz | fqzcomp5 > output.fqz5
fqzcomp5 -d input.fqz5 | gzip > output.fastq.gz
```

### Round-trip conversion
```bash
# Compress gzipped FASTQ
fqzcomp5 -5 data.fastq.gz data.fqz5

# Decompress to gzipped FASTQ
fqzcomp5 -d data.fqz5 data.out.fastq.gz

# Verify integrity
diff <(zcat data.fastq.gz) <(zcat data.out.fastq.gz)
```

## Implementation Details

### Key Components

1. **kseq.h**: Header-only FASTQ parser with gzip support
2. **load_seqs_kseq()**: New function to read FASTQ using kseq
3. **encode_gzip()**: Encoding function for gzipped input
4. **decode_gzip()** and **output_fastq_gzip()**: Decompression with gzipped output

### File Detection

The tool detects gzipped files by checking the `.gz` extension:
- Input files with `.gz` are read using `gzopen()`
- Output files with `.gz` are written using `gzwrite()`
- Plain files use standard `fopen()/fwrite()`

### Quality Score Handling

Quality scores are properly converted between internal format (0-based) and FASTQ format (Phred+33) in all code paths.

### Comment Preservation

FASTQ header comments (text after the first space) are preserved through compression/decompression cycles.

## Testing

Comprehensive tests verify:
- Plain FASTQ compression/decompression
- Gzipped input reading
- Gzipped output writing
- Round-trip conversions
- Quality score preservation
- Comment preservation
- Different compression levels (-1, -3, -5, -7, -9)

## Performance

The gzip support adds minimal overhead:
- Uses streaming I/O (no additional memory allocation)
- Leverages hardware-accelerated zlib when available
- Same compression ratios as plain FASTQ input
