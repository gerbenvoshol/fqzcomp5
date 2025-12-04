# Test Data for FQZComp5

This directory contains sample FASTQ files for testing the fqzcomp5 compression tool.

## Files

### sample.fastq
A small FASTQ file with 5 reads of varying sequences. This file is used to verify basic compression and decompression functionality.

## Running Tests

### Basic compression test:
```bash
# Compress the sample file
./fqzcomp5 test_data/sample.fastq test_data/sample.fqz5

# Decompress and verify
./fqzcomp5 -d test_data/sample.fqz5 test_data/sample.out.fastq
diff test_data/sample.fastq test_data/sample.out.fastq
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
