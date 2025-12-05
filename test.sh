#!/bin/bash

# Comprehensive test suite for fqzcomp5
# Tests roundtrip compression/decompression and read order preservation

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0

# Temporary directory for test outputs
TEST_DIR="/tmp/fqzcomp5_tests_$$"
mkdir -p "$TEST_DIR"

# Cleanup function
cleanup() {
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

# Function to print test results
print_result() {
    local test_name="$1"
    local status="$2"
    
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    
    if [ "$status" = "PASS" ]; then
        echo -e "${GREEN}[PASS]${NC} $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}[FAIL]${NC} $test_name"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

# Function to verify files are identical
verify_identical() {
    local file1="$1"
    local file2="$2"
    local test_name="$3"
    
    if diff -q "$file1" "$file2" > /dev/null 2>&1; then
        print_result "$test_name" "PASS"
        return 0
    else
        print_result "$test_name" "FAIL"
        echo "  Files differ: $file1 vs $file2"
        return 1
    fi
}

# Function to verify read order is preserved
verify_read_order() {
    local original="$1"
    local decompressed="$2"
    local test_name="$3"
    
    # Extract read names (lines starting with @)
    grep "^@" "$original" > "$TEST_DIR/orig_names.txt"
    grep "^@" "$decompressed" > "$TEST_DIR/decomp_names.txt"
    
    if diff -q "$TEST_DIR/orig_names.txt" "$TEST_DIR/decomp_names.txt" > /dev/null 2>&1; then
        print_result "$test_name" "PASS"
        return 0
    else
        print_result "$test_name" "FAIL"
        echo "  Read order differs between original and decompressed"
        return 1
    fi
}

echo "======================================"
echo "FQZComp5 Automated Test Suite"
echo "======================================"
echo ""

# Build the project if needed
if [ ! -f "./fqzcomp5" ]; then
    echo "Building fqzcomp5..."
    make -j4 > /dev/null 2>&1 || {
        echo -e "${RED}ERROR: Build failed${NC}"
        exit 1
    }
fi

echo "Running tests..."
echo ""

# Test 1: Basic single file roundtrip with default settings
echo "Test Group 1: Single File Compression"
echo "--------------------------------------"
./fqzcomp5 test_data/sample.fastq "$TEST_DIR/test1.fqz5" > /dev/null 2>&1
./fqzcomp5 -d "$TEST_DIR/test1.fqz5" "$TEST_DIR/test1.fastq" > /dev/null 2>&1
verify_identical test_data/sample.fastq "$TEST_DIR/test1.fastq" "Single file roundtrip (default)"
verify_read_order test_data/sample.fastq "$TEST_DIR/test1.fastq" "Single file read order preservation"

# Test 2-6: Different compression levels
for level in 1 3 5 7 9; do
    ./fqzcomp5 -$level test_data/sample.fastq "$TEST_DIR/test_level${level}.fqz5" > /dev/null 2>&1
    ./fqzcomp5 -d "$TEST_DIR/test_level${level}.fqz5" "$TEST_DIR/test_level${level}.fastq" > /dev/null 2>&1
    verify_identical test_data/sample.fastq "$TEST_DIR/test_level${level}.fastq" "Compression level -$level roundtrip"
done

echo ""
echo "Test Group 2: Paired-End Files"
echo "--------------------------------------"

# Test 7: Paired-end file compression with interleaving
./fqzcomp5 test_data/sample_R1.fastq test_data/sample_R2.fastq "$TEST_DIR/paired.fqz5" > /dev/null 2>&1
./fqzcomp5 -d "$TEST_DIR/paired.fqz5" "$TEST_DIR/paired_R1.fastq" "$TEST_DIR/paired_R2.fastq" > /dev/null 2>&1
verify_identical test_data/sample_R1.fastq "$TEST_DIR/paired_R1.fastq" "Paired-end R1 roundtrip"
verify_identical test_data/sample_R2.fastq "$TEST_DIR/paired_R2.fastq" "Paired-end R2 roundtrip"
verify_read_order test_data/sample_R1.fastq "$TEST_DIR/paired_R1.fastq" "Paired-end R1 read order preservation"
verify_read_order test_data/sample_R2.fastq "$TEST_DIR/paired_R2.fastq" "Paired-end R2 read order preservation"

# Test 8-12: Paired-end with different compression levels
for level in 1 3 5 7 9; do
    ./fqzcomp5 -$level test_data/sample_R1.fastq test_data/sample_R2.fastq "$TEST_DIR/paired_level${level}.fqz5" > /dev/null 2>&1
    ./fqzcomp5 -d "$TEST_DIR/paired_level${level}.fqz5" "$TEST_DIR/paired_R1_level${level}.fastq" "$TEST_DIR/paired_R2_level${level}.fastq" > /dev/null 2>&1
    verify_identical test_data/sample_R1.fastq "$TEST_DIR/paired_R1_level${level}.fastq" "Paired-end level -$level R1 roundtrip"
    verify_identical test_data/sample_R2.fastq "$TEST_DIR/paired_R2_level${level}.fastq" "Paired-end level -$level R2 roundtrip"
done

echo ""
echo "Test Group 3: Gzipped Input/Output"
echo "--------------------------------------"

# Create gzipped test files
gzip -c test_data/sample.fastq > "$TEST_DIR/sample.fastq.gz"
gzip -c test_data/sample_R1.fastq > "$TEST_DIR/sample_R1.fastq.gz"
gzip -c test_data/sample_R2.fastq > "$TEST_DIR/sample_R2.fastq.gz"

# Test 13: Gzipped input
./fqzcomp5 "$TEST_DIR/sample.fastq.gz" "$TEST_DIR/gzip_in.fqz5" > /dev/null 2>&1
./fqzcomp5 -d "$TEST_DIR/gzip_in.fqz5" "$TEST_DIR/gzip_in.fastq" > /dev/null 2>&1
verify_identical test_data/sample.fastq "$TEST_DIR/gzip_in.fastq" "Gzipped input roundtrip"

# Test 14: Gzipped output
./fqzcomp5 test_data/sample.fastq "$TEST_DIR/gzip_out.fqz5" > /dev/null 2>&1
./fqzcomp5 -d "$TEST_DIR/gzip_out.fqz5" "$TEST_DIR/gzip_out.fastq.gz" > /dev/null 2>&1
gunzip "$TEST_DIR/gzip_out.fastq.gz"
verify_identical test_data/sample.fastq "$TEST_DIR/gzip_out.fastq" "Gzipped output roundtrip"

# Test 15: Both gzipped input and output
./fqzcomp5 "$TEST_DIR/sample.fastq.gz" "$TEST_DIR/gzip_both.fqz5" > /dev/null 2>&1
./fqzcomp5 -d "$TEST_DIR/gzip_both.fqz5" "$TEST_DIR/gzip_both.fastq.gz" > /dev/null 2>&1
gunzip "$TEST_DIR/gzip_both.fastq.gz"
verify_identical test_data/sample.fastq "$TEST_DIR/gzip_both.fastq" "Gzipped input and output roundtrip"

# Test 16: Paired-end gzipped files
./fqzcomp5 "$TEST_DIR/sample_R1.fastq.gz" "$TEST_DIR/sample_R2.fastq.gz" "$TEST_DIR/paired_gzip.fqz5" > /dev/null 2>&1
./fqzcomp5 -d "$TEST_DIR/paired_gzip.fqz5" "$TEST_DIR/paired_gzip_R1.fastq.gz" "$TEST_DIR/paired_gzip_R2.fastq.gz" > /dev/null 2>&1
gunzip "$TEST_DIR/paired_gzip_R1.fastq.gz" "$TEST_DIR/paired_gzip_R2.fastq.gz"
verify_identical test_data/sample_R1.fastq "$TEST_DIR/paired_gzip_R1.fastq" "Paired-end gzipped R1 roundtrip"
verify_identical test_data/sample_R2.fastq "$TEST_DIR/paired_gzip_R2.fastq" "Paired-end gzipped R2 roundtrip"

echo ""
echo "Test Group 4: Special Options"
echo "--------------------------------------"

# Test 17: -p flag (output name on third line)
./fqzcomp5 test_data/sample.fastq "$TEST_DIR/plusname.fqz5" > /dev/null 2>&1
./fqzcomp5 -d -p "$TEST_DIR/plusname.fqz5" "$TEST_DIR/plusname.fastq" > /dev/null 2>&1

# Verify content is the same (excluding third line differences)
# Extract sequences and qualities only
grep -v "^[@+]" test_data/sample.fastq > "$TEST_DIR/orig_seqq.txt"
grep -v "^[@+]" "$TEST_DIR/plusname.fastq" > "$TEST_DIR/decomp_seqq.txt"
verify_identical "$TEST_DIR/orig_seqq.txt" "$TEST_DIR/decomp_seqq.txt" "Plus name flag (-p) sequences and qualities"

# Test 18: Multi-threaded compression
./fqzcomp5 -t 2 test_data/sample.fastq "$TEST_DIR/multithread.fqz5" > /dev/null 2>&1
./fqzcomp5 -d -t 2 "$TEST_DIR/multithread.fqz5" "$TEST_DIR/multithread.fastq" > /dev/null 2>&1
verify_identical test_data/sample.fastq "$TEST_DIR/multithread.fastq" "Multi-threaded (-t 2) roundtrip"

# Test 19: Custom block size
./fqzcomp5 -b 1K test_data/sample.fastq "$TEST_DIR/blocksize.fqz5" > /dev/null 2>&1
./fqzcomp5 -d "$TEST_DIR/blocksize.fqz5" "$TEST_DIR/blocksize.fastq" > /dev/null 2>&1
verify_identical test_data/sample.fastq "$TEST_DIR/blocksize.fastq" "Custom block size (-b 1K) roundtrip"

echo ""
echo "Test Group 5: File Format Validation"
echo "--------------------------------------"

# Test 20: Verify FQZ5 magic header
./fqzcomp5 test_data/sample.fastq "$TEST_DIR/header_check.fqz5" > /dev/null 2>&1
MAGIC=$(head -c 4 "$TEST_DIR/header_check.fqz5")
if [ "$MAGIC" = "FQZ5" ]; then
    print_result "FQZ5 magic header present" "PASS"
else
    print_result "FQZ5 magic header present" "FAIL"
fi

# Test 21: Verify index is present
if strings "$TEST_DIR/header_check.fqz5" | grep -q "FQZ5IDX"; then
    print_result "FQZ5 index present" "PASS"
else
    print_result "FQZ5 index present" "FAIL"
fi

echo ""
echo "Test Group 6: Integrity Checking"
echo "--------------------------------------"

# Test 22: Verify integrity check passes on good file
./fqzcomp5 test_data/sample.fastq "$TEST_DIR/integrity_test.fqz5" > /dev/null 2>&1
if ./fqzcomp5 --check "$TEST_DIR/integrity_test.fqz5" > /dev/null 2>&1; then
    print_result "Integrity check on valid file" "PASS"
else
    print_result "Integrity check on valid file" "FAIL"
fi

# Test 23: Verify integrity check detects corruption
cp "$TEST_DIR/integrity_test.fqz5" "$TEST_DIR/corrupted.fqz5"
# Corrupt the file by overwriting some bytes in the middle
dd if=/dev/zero of="$TEST_DIR/corrupted.fqz5" bs=1 count=10 seek=100 conv=notrunc > /dev/null 2>&1
if ! ./fqzcomp5 --check "$TEST_DIR/corrupted.fqz5" > /dev/null 2>&1; then
    print_result "Integrity check detects corruption" "PASS"
else
    print_result "Integrity check detects corruption" "FAIL"
fi

# Test 24: Verify --check with verbose mode
./fqzcomp5 test_data/sample.fastq "$TEST_DIR/verbose_check.fqz5" > /dev/null 2>&1
if ./fqzcomp5 --check -v "$TEST_DIR/verbose_check.fqz5" 2>&1 | grep -q "CRC OK"; then
    print_result "Integrity check verbose mode" "PASS"
else
    print_result "Integrity check verbose mode" "FAIL"
fi

# Test 25: Verify old format files report no CRC
if ./fqzcomp5 --check test_data/sample.fqz5 2>&1 | grep -q "no CRC"; then
    print_result "Old format files handled correctly" "PASS"
else
    print_result "Old format files handled correctly" "FAIL"
fi

echo ""
echo "======================================"
echo "Test Summary"
echo "======================================"
echo "Total tests: $TESTS_TOTAL"
echo -e "${GREEN}Passed: $TESTS_PASSED${NC}"
if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "${RED}Failed: $TESTS_FAILED${NC}"
    exit 1
else
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
fi
