# Test Data for FQZComp5

This directory contains sample FASTQ files for testing the fqzcomp5 compression tool.

## Files

### Single-End Files

#### sample.fastq
A small FASTQ file with 5 reads of varying sequences WITHOUT /1 or /2 suffixes. 
This file is used to verify basic single-end compression and decompression functionality.

#### single_with_suffix.fastq
A small FASTQ file with 5 reads WITH /1 suffixes. 
This tests the edge case of single-end files that have paired-end style naming.

### Paired-End Files

#### sample_R1.fastq and sample_R2.fastq
Paired-end FASTQ files (3 read pairs) WITH /1 and /2 suffixes.
These files test standard paired-end compression and decompression with explicit suffixes.

#### paired_R1_nosuffix.fastq and paired_R2_nosuffix.fastq
Paired-end FASTQ files (3 read pairs) WITHOUT /1 or /2 suffixes.
These files test the critical case where paired-end reads lack explicit suffixes - 
the compressor should add /1 and /2 flags during compression and restore them during decompression.
This was the bug fixed in PR #XX.

## Running Tests

### Basic single-end compression test:
```bash
# Compress the sample file (no suffixes)
./fqzcomp5 test_data/sample.fastq test_data/sample.fqz5

# Decompress and verify
./fqzcomp5 -d test_data/sample.fqz5 test_data/sample.out.fastq
diff test_data/sample.fastq test_data/sample.out.fastq
```

### Single-end with /1 suffix test:
```bash
# Test single-end file that has /1 suffix
./fqzcomp5 test_data/single_with_suffix.fastq test_data/single_suffix.fqz5
./fqzcomp5 -d test_data/single_suffix.fqz5 test_data/single_suffix.out.fastq
diff test_data/single_with_suffix.fastq test_data/single_suffix.out.fastq
```

### Paired-end compression test (with suffixes):
```bash
# Compress paired files with /1 /2 suffixes
./fqzcomp5 test_data/sample_R1.fastq test_data/sample_R2.fastq test_data/paired.fqz5

# Decompress and verify
./fqzcomp5 -d test_data/paired.fqz5 test_data/out_R1.fastq test_data/out_R2.fastq
diff test_data/sample_R1.fastq test_data/out_R1.fastq && echo "R1: OK"
diff test_data/sample_R2.fastq test_data/out_R2.fastq && echo "R2: OK"
```

### Paired-end compression test (WITHOUT suffixes - the critical bug case):
```bash
# Compress paired files without /1 /2 suffixes
# The compressor should add /1 and /2 during compression
./fqzcomp5 test_data/paired_R1_nosuffix.fastq test_data/paired_R2_nosuffix.fastq test_data/paired_nosuffix.fqz5

# Decompress - output files will have /1 and /2 added
./fqzcomp5 -d test_data/paired_nosuffix.fqz5 test_data/out_nosuffix_R1.fastq test_data/out_nosuffix_R2.fastq

# Verify the decompressed files have /1 and /2 suffixes added
head -1 test_data/out_nosuffix_R1.fastq | grep -q "/1" && echo "R1 suffix added: OK"
head -1 test_data/out_nosuffix_R2.fastq | grep -q "/2" && echo "R2 suffix added: OK"
```

### Test with different compression levels:
```bash
for level in 1 3 5 7 9; do
    echo "Testing level -$level"
    ./fqzcomp5 -$level test_data/sample.fastq test_data/sample_l$level.fqz5
    ./fqzcomp5 -d test_data/sample_l$level.fqz5 test_data/sample_l$level.out.fastq
    diff test_data/sample.fastq test_data/sample_l$level.out.fastq && echo "Level $level: OK"
done
```

### Test random access (using index):
The new file format includes an index that enables random access to blocks. You can verify the index is present:
```bash
# Check for index magic at end of file
tail -c 100 test_data/sample.fqz5 | hexdump -C | grep "FQZ5IDX"
```
