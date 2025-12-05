#!/bin/bash
# Test script to verify proper error handling for corrupted files

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Build the tool if needed
if [ ! -f "./fqzcomp5" ]; then
    make
fi

# Create a test file
echo "Creating test file..."
./fqzcomp5 -3 test_data/sample.fastq /tmp/test_error.fqz5 > /dev/null 2>&1

# Corrupt the file
echo "Corrupting file..."
dd if=/dev/zero of=/tmp/test_error.fqz5 bs=1 count=10 seek=50 conv=notrunc 2>&1 | grep -v records

# Test that decompression fails gracefully without segfault
echo "Testing error handling..."
if ./fqzcomp5 -d /tmp/test_error.fqz5 /tmp/test_error_out.fastq 2>&1 | grep -q "ERROR:"; then
    echo -e "${GREEN}PASS: Error handling works correctly (no segfault)${NC}"
    rm -f /tmp/test_error.fqz5 /tmp/test_error_out.fastq
    exit 0
else
    echo -e "${RED}FAIL: No error message detected${NC}"
    rm -f /tmp/test_error.fqz5 /tmp/test_error_out.fastq
    exit 1
fi
