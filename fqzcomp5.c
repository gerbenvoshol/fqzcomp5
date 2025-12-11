/* Tests for fqz codec */
/*
 * Copyright (c) 2019,2020,2022 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
File format (Version 1.1):

[Header]
8    Magic number: "FQZ5\001\001\000\000"  (FQZ5 followed by version 1.1.0)
8    Index offset (0 if no index)

[Block]*  Zero or more blocks of records

4    Block size (total bytes in this block, excluding this 4-byte field)
4    Num records
4    Block CRC32 checksum (CRC of all data from num_records onwards)

1    Name strategy
4    NU: uncompressed name length
4    NC: compressed name length
NC   Name compressed data

1    Read length strategy
?    Read length data

1    Sequence strategy (bits 0..2), both_strands (bit 3), level (4..7)
4    SU: Uncompressed sequence size
4    SC: Compressed sequence size
SC   Compressed sequence data

1    Quality strategy
4    QU: Uncompressed quality size
4    QC: Compressed quality size
QC   Compressed quality data

[Index] (optional, at end of file)
8    Magic: "FQZ5IDX\000"
4    Number of blocks
For each block:
  8  File offset of block
  4  Uncompressed size (total bases)
  4  Number of records in block
4    Index CRC32 checksum

[Trailer] (optional, at end of file after index)
8    Magic: "FQZ5END\000"
4    Overall file CRC32 (CRC of all block data)
4    Number of blocks

Version 1.0 files (without CRC) are still supported for backward compatibility.

 */

// TODO
// - Split aux tags into own data series using CRAM's TL + per tag.

// - Seq encoding using STR + copy-number?
//   Couldn't get this to work well though, even though it sounds like an
//   ideal thing for ONT/PB-CLR.  Maybe let it kick in after so many bases
//   in a homopolymer?

// - Removal of lengths from fqzqual stream

// - Entropy encoding of read length stream

// - Reuse of memory buffers for speed

// - Also check why memory usage is so high.  Is it all in data models?

// - Increase sizes of fqzcomp contexts?  Not so useful for CRAM, but maybe
//   it's still beneficial to go beyond 16-bit.

// - Improve fqzcomp tables to permit any mapping rather than monotonic.

// - Distinguish explicit method opts (-s1 -S13B -q1 -Q2 etc) from auto
//   picked options (-3, -5) which use auto-selected metrics

// - We don't need large input fastq blocks to do large block compression.
//   We could separate model reset from block boundaries, so blocks are small
//   but reset boundaries less common.

// - Flag for "+" (3rd record) duplicating first "@".

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>
#include <pthread.h>
#include <sys/time.h>
#include <zlib.h>

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#include "htscodecs/varint.h"
#include "htscodecs/fqzcomp_qual.h"
#include "htscodecs/tokenise_name3.h"
#include "htscodecs/rANS_static4x16.h"
#include "htscodecs/varint.h"
#include "thread_pool.h"
#include "lzp16e.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#define BLK_SIZE 512*1000000

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

// Review metrics stats every X blocks for Y trials
#define METRICS_REVIEW 100
#define METRICS_TRIAL 3

// File format constants
#define FQZ5_MAGIC_V10 "FQZ5\001\000\000\000"  // Magic + version 1.0.0 (old)
#define FQZ5_MAGIC "FQZ5\001\001\000\000"      // Magic + version 1.1.0 (current)
#define FQZ5_MAGIC_LEN 8
#define FQZ5_INDEX_MAGIC "FQZ5IDX\000"
#define FQZ5_INDEX_MAGIC_LEN 8
#define FQZ5_TRAILER_MAGIC "FQZ5END\000"
#define FQZ5_TRAILER_MAGIC_LEN 8

// Index entry for each block
typedef struct {
    uint64_t offset;       // File offset of block start
    uint32_t usize;        // Uncompressed size (total bases)
    uint32_t nrecords;     // Number of records in block
} index_entry;

// File index structure
typedef struct {
    uint32_t nblocks;      // Number of blocks
    index_entry *entries;  // Array of index entries
} fqz5_index;

typedef enum {
    SEC_NAME,
    SEC_LEN,
    SEC_SEQ,
    SEC_QUAL,
    SEC_LAST
} sections;

// TODO: add method costs too
typedef enum {
    // general purpose
    RANS0=1, RANS1, RANS64, RANS65, RANS128, RANS129, RANS192, RANS193,
    RANSXN1,

    // LZP; differing min lengths?  Make len part of format?
    LZP3,
    TLZP3, // name lzp; fold into LZP3 with better format structuring.
    // TODO LZP2, LZP4, LZP16? Needs storing in byte stream too

    // Name specific; may just use arg.slevel and ignore multiplicity here?
    // Do we ever want to gather starts on multiple compression levels to
    // judge if worth it?  Maybe, but only for highest probably.
    TOK3_3, TOK3_5, TOK3_7, TOK3_9,
    TOK3_3_LZP, TOK3_5_LZP, TOK3_7_LZP, TOK3_9_LZP,

    // Seq
    SEQ10, SEQ12, SEQ12B, SEQ13B, SEQ14B, SEQ_CUSTOM,

    // Qual
    FQZ0, FQZ1, FQZ2, FQZ3, FQZ4,

    M_LAST,
} methods;

int method_costs[] = {
    1.0, 1.1, 1.1, 1.2, 0.5, 0.7, 1.0, 1.2, 1.2, // various RANS

    1.3, // LZP3
    1.3, // TLZP3

    1.3, 1.4, 1.5, 1.6, // TOK3
    1.3, 1.4, 1.5, 1.6, // TOK3+LZP

    1.7, 1.8, 1.9, 2.0, 2.2, 1.5, // FQZ-SEQ

    1.3, 1.3, 1.3, 1.3, 1.3, // FQZ-QUAL
};

typedef struct {
    uint64_t usize[M_LAST], csize[M_LAST];   // current accumulated sizes
    int review, trial;
    int count[M_LAST];
} metrics;

pthread_mutex_t metric_m = PTHREAD_MUTEX_INITIALIZER;
static uint32_t method_avail[SEC_LAST];
static methods method_used[SEC_LAST];
static metrics stats[SEC_LAST];

typedef struct {
    int num_records;                 // number of fastq entries
    char *name_buf;                  // concatenated names, \0 separator
    char *seq_buf;                   // concatenated seq,   no separator
    char *qual_buf;                  // concatenated qual,  no separator
    int *name;                       // index into name_buf
    int *seq;                        // index into seq_buf
    int *qual;                       // index into qual_buf
    unsigned int *len;               // sequence length
    unsigned int *flag;              // READ1/READ2 parsed from name
    int name_len, seq_len, qual_len; // used size of _buf above
    int name_sz,  seq_sz,  qual_sz;  // alloced size of _buf above
    int fixed_len;                   // length of each seq, 0 if not fixed
    int is_fasta;                    // 1 if FASTA format (no quality), 0 if FASTQ
} fastq;

fastq *fastq_alloc(int nr) {
    fastq *fq = calloc(1, sizeof(*fq));

    fq->num_records = nr;
    fq->name = calloc(nr, sizeof(char *));
    fq->seq  = calloc(nr, sizeof(char *));
    fq->qual = calloc(nr, sizeof(char *));
    fq->len  = calloc(nr, sizeof(int));
    fq->flag = calloc(nr, sizeof(int));

    return fq;
}

void fastq_free(fastq *fq) {
    if (!fq)
	return;
    free(fq->name_buf);
    free(fq->seq_buf);
    free(fq->qual_buf);
    free(fq->name);
    free(fq->seq);
    free(fq->qual);
    free(fq->len);
    free(fq->flag);
    free(fq);
}

#define goto if (fprintf(stderr, "ERR %s:%d\n", __FILE__, __LINE__)) goto
fastq *load_seqs(char *in, int blk_size, int *last_offset) {
    fastq *fq = calloc(1, sizeof(*fq));
    if (!fq)
	goto err;
    size_t name_sz = blk_size/10;
    size_t seq_sz  = blk_size/2;
    size_t qual_sz = blk_size/2;
    char *name_buf = fq->name_buf = malloc(name_sz);
    char *seq_buf  = fq->seq_buf  = malloc(seq_sz);
    char *qual_buf = fq->qual_buf = malloc(qual_sz);
    if (!name_buf || !seq_buf || !qual_buf)
	goto err;
    int c, i = 0, nr = 0, ar = 0;
    int last_name = -1;
    fq->fixed_len = -1;

    int name_i = 0, seq_i = 0, qual_i = 0;
    int last_start = 0;
    while (i < blk_size) {
	if (nr >= ar) {
	    ar = ar ? (ar << 1) : 10000;  // 2x growth, more cache-friendly
	    fq->name = realloc(fq->name, ar*sizeof(char *));
	    fq->seq  = realloc(fq->seq , ar*sizeof(char *));
	    fq->qual = realloc(fq->qual, ar*sizeof(char *));
	    fq->len  = realloc(fq->len,  ar*sizeof(int));
	    fq->flag = realloc(fq->flag, ar*sizeof(int));
	}

	// @name
	fq->name[nr] = name_i;
	c = in[i++];
	if (c != '@')
	    goto err;

	int name_i_ = name_i; // tmp copy so we can unwind a partial decode
	while (i < blk_size && (c = in[i++]) && c != '\n') {
	    if (name_i_+1 >= name_sz) {
		name_sz = name_sz ? (name_sz << 1) : 1000;  // 2x growth
		name_buf = fq->name_buf = realloc(fq->name_buf, name_sz);
	    }
	    name_buf[name_i_++] = c;
	}
	if (i == blk_size)
	    break;

	name_buf[name_i_++] = 0;

	int flag = 0;
	if (name_i_ > 3 &&
	    name_buf[name_i_-2] == '2' &&
	    name_buf[name_i_-3] == '/')
	    flag = FQZ_FREAD2;
	if (last_name >= 0 &&
	    strcmp(fq->name_buf + fq->name[nr], fq->name_buf + last_name) == 0)
	    flag = FQZ_FREAD2;
	fq->flag[nr] = flag;
	last_name = fq->name[nr];

	// seq
	fq->seq[nr] = seq_i;
	int len = seq_i;
	int seq_i_ = seq_i;
	while (i < blk_size && (c = in[i++]) && c != '\n') {
	    if (seq_i_ >= seq_sz) {
		// very unlikely given blk_size/2 starting point,
		// but not impossible.
		seq_sz = seq_sz ? (seq_sz << 1) : 1000;  // 2x growth
		seq_buf = fq->seq_buf = realloc(fq->seq_buf, seq_sz);
	    }
	    seq_buf[seq_i_++] = c;
	}
	if (i == blk_size)
	    break;
	fq->len[nr] = seq_i_ - len;
	//seq_buf[seq_i_++] = '\n'; // use instead of length terminator?

	if (fq->fixed_len == -1)
	    fq->fixed_len = fq->len[nr];
	else if (fq->fixed_len > 0)
	    if (fq->fixed_len != fq->len[nr])
		fq->fixed_len = 0;

	// +(name)
	if (i < blk_size && (c = in[i++]) != '+')
	    goto err;
	while (i < blk_size && (c = in[i++]) && c != '\n')
	    ;
	if (i == blk_size)
	    break;

	// qual
	fq->qual[nr] = qual_i;
	len = qual_i;
	int qual_i_ = qual_i;
	while (i < blk_size && (c = in[i++]) && c != '\n') {
	    if (qual_i_ >= qual_sz) {
		// very unlikely given blk_size/2 starting point
		qual_sz = qual_sz ? (qual_sz << 1) : 1000;  // 2x growth
		qual_buf = fq->qual_buf = realloc(fq->qual_buf, qual_sz);
	    }
	    qual_buf[qual_i_++] = c-33;
	}

	if (fq->len[nr] != qual_i_ - len) {
	    if (i == blk_size)
		break;

	    goto err;
	} else if (i == blk_size && c != '\n')
	    break;

	name_i = name_i_;
	seq_i  = seq_i_;
	qual_i = qual_i_;

	last_start = i;
	nr++;
    }

    // FIXME: need to deal with case where block size is smaller than
    // a single record!  For now that puts a limit on smallest block size.
    *last_offset = last_start; // so we can continue for next block
    
    fq->name_len = name_i;
    fq->seq_len  = seq_i;
    fq->qual_len = qual_i;
    fq->num_records = nr;

    // Reduce memory wastage
    fq->name_buf = realloc(fq->name_buf, name_i);
    fq->seq_buf  = realloc(fq->seq_buf,  seq_i);
    fq->qual_buf = realloc(fq->qual_buf, qual_i);

    return fq;

 err:
    fprintf(stderr, "Failed to load fastq input\n");
    fastq_free(fq);

    return NULL;
}
#undef goto

// Load sequences using kseq.h (supports gzipped files)
fastq *load_seqs_kseq(gzFile fp, int blk_size, int *eof_flag) {
    // Static variables to maintain state across calls for the same file
    static kseq_t *seq = NULL;
    static gzFile last_fp = NULL;
    static int have_buffered = 0;
    
    // Reinitialize if file handle changed or first call
    if (seq == NULL || fp != last_fp) {
        if (seq != NULL) {
            kseq_destroy(seq);
        }
        seq = kseq_init(fp);
        last_fp = fp;
        have_buffered = 0;
    }
    
    fastq *fq = calloc(1, sizeof(*fq));
    if (!fq)
        goto err;
    
    size_t name_sz = blk_size/10;
    size_t seq_sz  = blk_size/2;
    size_t qual_sz = blk_size/2;
    char *name_buf = fq->name_buf = malloc(name_sz);
    char *seq_buf  = fq->seq_buf  = malloc(seq_sz);
    char *qual_buf = fq->qual_buf = malloc(qual_sz);
    if (!name_buf || !seq_buf || !qual_buf)
        goto err;
    
    int nr = 0, ar = 0;
    int last_name = -1;
    fq->fixed_len = -1;
    
    int name_i = 0, seq_i = 0, qual_i = 0;
    int total_size = 0;
    
    int l = 0;
    
    // Process buffered record from previous block if we have one
    if (have_buffered) {
        have_buffered = 0;
        l = 0;  // Signal we have a valid record
        goto process_record;
    }
    
    while ((l = kseq_read(seq)) >= 0) {
    process_record:
        // Check if we would exceed block size
        int record_size = seq->name.l + 1 + seq->seq.l + seq->qual.l;
        if (total_size > 0 && total_size + record_size > blk_size) {
            // Block is full - buffer this record for next block
            have_buffered = 1;
            break;
        }
        
        total_size += record_size;
        
        if (nr >= ar) {
            ar = ar ? (ar << 1) : 10000;  // 2x growth, more cache-friendly
            fq->name = realloc(fq->name, ar*sizeof(char *));
            fq->seq  = realloc(fq->seq , ar*sizeof(char *));
            fq->qual = realloc(fq->qual, ar*sizeof(char *));
            fq->len  = realloc(fq->len,  ar*sizeof(int));
            fq->flag = realloc(fq->flag, ar*sizeof(int));
        }
        
        // Store name (and comment if present)
        fq->name[nr] = name_i;
        int total_name_len = seq->name.l;
        if (seq->comment.l > 0) {
            total_name_len += 1 + seq->comment.l; // space + comment
        }
        
        if (name_i + total_name_len + 1 >= name_sz) {
            // Calculate exact required size, ensuring at least 2x growth
            size_t required = name_i + total_name_len + 1;
            size_t new_sz = name_sz ? (name_sz << 1) : 1000;
            if (new_sz < required)
                new_sz = required;
            name_sz = new_sz;
            name_buf = fq->name_buf = realloc(fq->name_buf, name_sz);
        }
        memcpy(name_buf + name_i, seq->name.s, seq->name.l);
        name_i += seq->name.l;
        
        // Add comment if present
        if (seq->comment.l > 0) {
            name_buf[name_i++] = ' ';
            memcpy(name_buf + name_i, seq->comment.s, seq->comment.l);
            name_i += seq->comment.l;
        }
        
        name_buf[name_i++] = 0;
        
        int flag = 0;
        if (seq->name.l > 1 &&
            name_buf[name_i-2] == '2' &&
            name_buf[name_i-3] == '/')
            flag = FQZ_FREAD2;
        if (last_name >= 0 &&
            strcmp(fq->name_buf + fq->name[nr], fq->name_buf + last_name) == 0)
            flag = FQZ_FREAD2;
        fq->flag[nr] = flag;
        last_name = fq->name[nr];
        
        // Store sequence
        fq->seq[nr] = seq_i;
        if (seq_i + seq->seq.l >= seq_sz) {
            // Calculate exact required size, ensuring at least 2x growth
            size_t required = seq_i + seq->seq.l;
            size_t new_sz = seq_sz ? (seq_sz << 1) : 1000;
            if (new_sz < required)
                new_sz = required;
            seq_sz = new_sz;
            seq_buf = fq->seq_buf = realloc(fq->seq_buf, seq_sz);
        }
        memcpy(seq_buf + seq_i, seq->seq.s, seq->seq.l);
        seq_i += seq->seq.l;
        fq->len[nr] = seq->seq.l;
        
        if (fq->fixed_len == -1)
            fq->fixed_len = fq->len[nr];
        else if (fq->fixed_len > 0)
            if (fq->fixed_len != fq->len[nr])
                fq->fixed_len = 0;
        
        // Store quality (only for FASTQ)
        fq->qual[nr] = qual_i;
        if (seq->qual.l > 0) {
            // FASTQ format
            if (qual_i + seq->qual.l >= qual_sz) {
                // Calculate exact required size, ensuring at least 2x growth
                size_t required = qual_i + seq->qual.l;
                size_t new_sz = qual_sz ? (qual_sz << 1) : 1000;
                if (new_sz < required)
                    new_sz = required;
                qual_sz = new_sz;
                qual_buf = fq->qual_buf = realloc(fq->qual_buf, qual_sz);
            }
            for (int i = 0; i < seq->qual.l; i++)
                qual_buf[qual_i++] = seq->qual.s[i] - 33;
        }
        
        if (seq->qual.l > 0) {
            // FASTQ format
            if (fq->len[nr] != seq->qual.l) {
                fprintf(stderr, "Sequence and quality length mismatch\n");
                kseq_destroy(seq);
                goto err;
            }
        } else {
            // FASTA format
            if (nr == 0) {
                fq->is_fasta = 1;
            }
        }
        
        nr++;
    }
    
    // Check EOF status
    if (l == -1) {
        *eof_flag = 1;  // Normal EOF
        // Clean up static kseq on EOF
        kseq_destroy(seq);
        seq = NULL;
        have_buffered = 0;
    } else if (l < -1) {
        fprintf(stderr, "Error reading FASTQ file (code %d)\n", l);
        // Clean up on error
        kseq_destroy(seq);
        seq = NULL;
        have_buffered = 0;
        goto err;
    }
    
    fq->name_len = name_i;
    fq->seq_len  = seq_i;
    fq->qual_len = qual_i;
    fq->num_records = nr;
    
    // Reduce memory wastage
    fq->name_buf = realloc(fq->name_buf, name_i > 0 ? name_i : 1);
    fq->seq_buf  = realloc(fq->seq_buf,  seq_i > 0 ? seq_i : 1);
    fq->qual_buf = realloc(fq->qual_buf, qual_i > 0 ? qual_i : 1);
    
    return fq;

 err:
    fprintf(stderr, "Failed to load fastq input with kseq\n");
    fastq_free(fq);
    // Clean up static state on error to prevent memory leaks and invalid state
    if (seq != NULL) {
        kseq_destroy(seq);
        seq = NULL;
        last_fp = NULL;
        have_buffered = 0;
    }
    return NULL;
}

// Load interleaved sequences from two paired FASTQ files
// Reads alternately from fp1 and fp2 to create an interleaved stream
fastq *load_seqs_interleaved(gzFile fp1, gzFile fp2, int blk_size, int *eof_flag) {
    // Static variables to maintain state across calls for the same files
    static kseq_t *seq1 = NULL;
    static kseq_t *seq2 = NULL;
    static gzFile last_fp1 = NULL;
    static gzFile last_fp2 = NULL;
    static int have_buffered = 0;
    
    // Reinitialize if file handles changed or first call
    if (seq1 == NULL || seq2 == NULL || fp1 != last_fp1 || fp2 != last_fp2) {
        if (seq1 != NULL) {
            kseq_destroy(seq1);
        }
        if (seq2 != NULL) {
            kseq_destroy(seq2);
        }
        seq1 = kseq_init(fp1);
        seq2 = kseq_init(fp2);
        last_fp1 = fp1;
        last_fp2 = fp2;
        have_buffered = 0;
    }
    
    fastq *fq = calloc(1, sizeof(*fq));
    if (!fq)
        goto err;
    
    size_t name_sz = blk_size/10;
    size_t seq_sz  = blk_size/2;
    size_t qual_sz = blk_size/2;
    char *name_buf = fq->name_buf = malloc(name_sz);
    char *seq_buf  = fq->seq_buf  = malloc(seq_sz);
    char *qual_buf = fq->qual_buf = malloc(qual_sz);
    if (!name_buf || !seq_buf || !qual_buf)
        goto err;
    
    int nr = 0, ar = 0;
    int last_name = -1;
    fq->fixed_len = -1;
    
    int name_i = 0, seq_i = 0, qual_i = 0;
    int total_size = 0;
    
    int l1 = 0, l2 = 0;
    
    // Process buffered record pair from previous block if we have one
    if (have_buffered) {
        have_buffered = 0;
        l1 = 0;  // Signal we have valid records
        l2 = 0;
        goto process_pair;
    }
    
    // Read pairs alternately until block is full or EOF
    while (1) {
        // Read from file 1 (R1)
        l1 = kseq_read(seq1);
        if (l1 < 0)
            break;
        
        // Read from file 2 (R2)
        l2 = kseq_read(seq2);
        if (l2 < 0) {
            fprintf(stderr, "Unpaired read detected: R2 file ended before R1\n");
            // Clean up on error
            kseq_destroy(seq1);
            kseq_destroy(seq2);
            seq1 = NULL;
            seq2 = NULL;
            have_buffered = 0;
            goto err;
        }
        
    process_pair:
        // Check if we would exceed block size with both reads
        int record_size1 = seq1->name.l + 1 + seq1->seq.l + seq1->qual.l;
        int record_size2 = seq2->name.l + 1 + seq2->seq.l + seq2->qual.l;
        if (total_size > 0 && total_size + record_size1 + record_size2 > blk_size) {
            // Block is full - buffer this pair for next block
            have_buffered = 1;
            break;
        }
        
        total_size += record_size1 + record_size2;
        
        // Process both R1 and R2 reads
        for (int pair = 0; pair < 2; pair++) {
            kseq_t *seq = (pair == 0) ? seq1 : seq2;
            
            if (nr >= ar) {
                ar = ar*1.5 + 10000;
                int *new_name = realloc(fq->name, ar*sizeof(char *));
                int *new_seq  = realloc(fq->seq , ar*sizeof(char *));
                int *new_qual = realloc(fq->qual, ar*sizeof(char *));
                unsigned int *new_len  = realloc(fq->len,  ar*sizeof(int));
                unsigned int *new_flag = realloc(fq->flag, ar*sizeof(int));
                
                if (!new_name || !new_seq || !new_qual || !new_len || !new_flag) {
                    fprintf(stderr, "Failed to reallocate arrays for interleaved reads\n");
                    kseq_destroy(seq1);
                    kseq_destroy(seq2);
                    goto err;
                }
                
                fq->name = new_name;
                fq->seq  = new_seq;
                fq->qual = new_qual;
                fq->len  = new_len;
                fq->flag = new_flag;
            }
            
            // Store name (and comment if present)
            fq->name[nr] = name_i;
            int total_name_len = seq->name.l;
            if (seq->comment.l > 0) {
                total_name_len += 1 + seq->comment.l; // space + comment
            }
            
            if (name_i + total_name_len + 1 >= name_sz) {
                name_sz = name_sz * 1.5 + total_name_len + 1000;
                name_buf = fq->name_buf = realloc(fq->name_buf, name_sz);
            }
            memcpy(name_buf + name_i, seq->name.s, seq->name.l);
            name_i += seq->name.l;
            
            // Add comment if present
            if (seq->comment.l > 0) {
                name_buf[name_i++] = ' ';
                memcpy(name_buf + name_i, seq->comment.s, seq->comment.l);
                name_i += seq->comment.l;
            }
            
            name_buf[name_i++] = 0;
            
            // Set flag for R2 reads
            int flag = (pair == 1) ? FQZ_FREAD2 : 0;
            fq->flag[nr] = flag;
            if (pair == 0)
                last_name = fq->name[nr];
            
            // Store sequence
            fq->seq[nr] = seq_i;
            if (seq_i + seq->seq.l >= seq_sz) {
                seq_sz = seq_sz * 1.5 + seq->seq.l + 1000;
                seq_buf = fq->seq_buf = realloc(fq->seq_buf, seq_sz);
            }
            memcpy(seq_buf + seq_i, seq->seq.s, seq->seq.l);
            seq_i += seq->seq.l;
            fq->len[nr] = seq->seq.l;
            
            if (fq->fixed_len == -1)
                fq->fixed_len = fq->len[nr];
            else if (fq->fixed_len > 0)
                if (fq->fixed_len != fq->len[nr])
                    fq->fixed_len = 0;
            
            // Store quality (only for FASTQ)
            fq->qual[nr] = qual_i;
            if (seq->qual.l > 0) {
                // FASTQ format
                if (qual_i + seq->qual.l >= qual_sz) {
                    qual_sz = qual_sz * 1.5 + seq->qual.l + 1000;
                    qual_buf = fq->qual_buf = realloc(fq->qual_buf, qual_sz);
                }
                for (int i = 0; i < seq->qual.l; i++)
                    qual_buf[qual_i++] = seq->qual.s[i] - 33;
                
                if (fq->len[nr] != seq->qual.l) {
                    fprintf(stderr, "Sequence and quality length mismatch\n");
                    // Clean up on error
                    kseq_destroy(seq1);
                    kseq_destroy(seq2);
                    seq1 = NULL;
                    seq2 = NULL;
                    have_buffered = 0;
                    goto err;
                }
            } else {
                // FASTA format - no quality scores
                if (nr == 0) {
                    fq->is_fasta = 1;
                }
            }
            
            nr++;
        }
    }
    
    // Check EOF status and clean up if needed
    if (l1 == -1 && l2 == -1) {
        *eof_flag = 1;  // Normal EOF
        // Clean up static kseq on EOF
        kseq_destroy(seq1);
        kseq_destroy(seq2);
        seq1 = NULL;
        seq2 = NULL;
        have_buffered = 0;
    } else if (l1 < -1 || l2 < -1) {
        fprintf(stderr, "Error reading FASTQ files (codes %d, %d)\n", l1, l2);
        // Clean up on error
        kseq_destroy(seq1);
        kseq_destroy(seq2);
        seq1 = NULL;
        seq2 = NULL;
        have_buffered = 0;
        goto err;
    }
    // Note: If l1 >= 0 and l2 >= 0, we broke due to block size and will continue next time
    
    fq->name_len = name_i;
    fq->seq_len  = seq_i;
    fq->qual_len = qual_i;
    fq->num_records = nr;
    
    // Reduce memory wastage
    fq->name_buf = realloc(fq->name_buf, name_i > 0 ? name_i : 1);
    fq->seq_buf  = realloc(fq->seq_buf,  seq_i > 0 ? seq_i : 1);
    fq->qual_buf = realloc(fq->qual_buf, qual_i > 0 ? qual_i : 1);
    
    return fq;

 err:
    fprintf(stderr, "Failed to load interleaved fastq input\n");
    fastq_free(fq);
    // Clean up static state on error to prevent memory leaks and invalid state
    if (seq1 != NULL) {
        kseq_destroy(seq1);
        seq1 = NULL;
    }
    if (seq2 != NULL) {
        kseq_destroy(seq2);
        seq2 = NULL;
    }
    last_fp1 = NULL;
    last_fp2 = NULL;
    have_buffered = 0;
    return NULL;
}

#if 0
static uint64_t manual_strats[10] = {0};
static int manual_nstrat = 0;

/*
 * Manually specified strategies held in global manual_strats[].
 */
static inline
int fqz_manual_parameters(fqz_gparams *gp,
			  fqz_slice *s,
			  unsigned char *in,
			  size_t in_size) {
    int i, p;
    int dsqr[] = {
	0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
    };

    gp->vers = FQZ_VERS;
    gp->nparam = manual_nstrat;
    gp->gflags = GFLAG_MULTI_PARAM | GFLAG_HAVE_STAB;
    for (i = 0; i < 256; i++)
	gp->stab[i] = 0;

    // Fill these out later
    gp->max_sel = 0;
    gp->max_sym = 0;
    gp->p = malloc(gp->nparam * sizeof(*gp->p));

    for (p = 0; p < gp->nparam; p++) {
	fqz_param *pm = &gp->p[p];
	uint64_t st = manual_strats[p];

	pm->boff   = st & 15; st >>= 4;
	pm->bloc   = st & 15; st >>= 4;
	pm->bbits  = st & 15; st >>= 4;
	pm->do_qa  = st & 15; st >>= 4;
	pm->do_r2  = st & 15; st >>= 4;
	pm->dloc   = st & 15; st >>= 4;
	pm->ploc   = st & 15; st >>= 4;
	pm->sloc   = st & 15; st >>= 4;
	pm->qloc   = st & 15; st >>= 4;
	pm->dshift = st & 15; st >>= 4;
	pm->dbits  = st & 15; st >>= 4;
	pm->pshift = st & 15; st >>= 4;
	pm->pbits  = st & 15; st >>= 4;
	pm->qshift = st & 15; st >>= 4;
	pm->qbits  = st & 15; st >>= 4;

	// Gather some stats, as per qual_stats func.
	// r in rec count.
	// i = index to in[]
	// j = index within this rec
	uint32_t qhist[256] = {0};

	// qual stats for seqs using this parameter only
	fqz_qual_stats(s, in, in_size, pm, qhist, p);
	int max_sel = pm->max_sel;

	// Update max_sel running total. Eg with 4 sub-params:
	//
	// sel    param no.   => new
	// 0      0              0
	// 0/1    1              1,2
	// 0/1    2              3,4
	// 0      3              5
	for (i = gp->max_sel; i < gp->max_sel + max_sel+1; i++)
	    gp->stab[i] = p;
	gp->max_sel += max_sel+1;

	pm->fixed_len = pm->fixed_len > 0;
	pm->use_qtab = 0;  // unused by current encoder
	pm->store_qmap = pm->nsym <= 8;

	// Adjust parameters based on quality stats.
	// FIXME:  dup from fqz_pick_parameters.
	for (i = 0; i < sizeof(dsqr)/sizeof(*dsqr); i++)
	    if (dsqr[i] > (1<<pm->dbits)-1)
		dsqr[i] = (1<<pm->dbits)-1;

	if (pm->store_qmap) {
	    int j;
	    for (i = j = 0; i < 256; i++)
		if (qhist[i])
		    pm->qmap[i] = j++;
		else
		    pm->qmap[i] = INT_MAX;
	    pm->max_sym = pm->nsym;
	} else {
	    pm->nsym = 255;
	    for (i = 0; i < 256; i++)
		pm->qmap[i] = i;
	}
	if (gp->max_sym < pm->max_sym)
	    gp->max_sym = pm->max_sym;

	// Produce ptab from pshift.
	if (pm->qbits) {
	    for (i = 0; i < 256; i++) {
		pm->qtab[i] = i; // 1:1
		//pm->qtab[i] = (1<<pm->qshift)-i; // 1:1
		//pm->qtab[i] = i/4;
		//pm->qtab[i] = ((1<<pm->qshift)-i)/3;

		// Alternative mappings:
		//qtab[i] = i > 30 ? MIN(max_sym,i)-15 : i/2;  // eg for 9827 BAM
	    }
#if 0
	    // qtab for PacBio CCS data; saves 3%
	    for (i='~'-33; i<256; i++) {
		pm->qtab[i] = 24+i-('~'-33);
	    }

	    int x = 0;
	    for (i = 0; i < 1; i++)
		pm->qtab[i] = x,x++;
	    for (;i < '~'-33; i++)
		//pm->qtab[i] = x,x+=(i%4==0);
		pm->qtab[i] = x,x+=(i%4==0);
	    x++;
	    for (;i <= '~'-33; i++)
		pm->qtab[i] = x,x++;
	    for (;i < 256; i++)
		pm->qtab[i] = x,x+=(i%4==0);
#endif

//	    for (i = 0; i < 128; i++) {
//		for (x = i+1; x < 128; x++) {
//		    if (pm->qtab[i] != pm->qtab[x])
//			break;
//		}
//		x--;
//		if (i==x)
//		    fprintf(stderr, "%d:%d ", pm->qtab[i], i);
//		else {
//		    fprintf(stderr, "%d:%d-%d ", pm->qtab[i], i, x);
//		    i=x;
//		}
//	    }
//	    fprintf(stderr, "\n");

	    //pm->qtab['~'-33]=32;

	    // pm->use_qtab = 1;
//	    for (i = 0; i <= 2 ; i++) pm->qtab[i] = 0;
//	    for (     ; i <= 12; i++) pm->qtab[i] = 1;
//	    for (     ; i <= 18; i++) pm->qtab[i] = 2;
//	    for (     ; i <= 36; i++) pm->qtab[i] = 3;
	}
	//pm->use_qtab = 0;
	pm->qmask = (1<<pm->qbits)-1;

	if (pm->pbits) {
	    for (i = 0; i < 1024; i++)
		pm->ptab[i] = MIN((1<<pm->pbits)-1, i>>pm->pshift);

//	    for (i = 0; i < 1024; i++)
//		pm->ptab[i] = MIN((1<<pm->pbits)-1, i < 10 ? i : 10 + i/3);

	    // Alternatively via analysis of quality distributions we
	    // may select a bunch of positions that are special and
	    // have a non-uniform ptab[].
	    // Manual experimentation on a NovaSeq run saved 2.8% here.
	}

	if (pm->dbits) {
	    for (i = 0; i < 256; i++)
		pm->dtab[i] = dsqr[MIN(sizeof(dsqr)/sizeof(*dsqr)-1, i>>pm->dshift)];
	}

	pm->use_ptab = (pm->pbits > 0);
	pm->use_dtab = (pm->dbits > 0);

	pm->pflags =
	    (pm->use_qtab   ?PFLAG_HAVE_QTAB :0)|
	    (pm->use_dtab   ?PFLAG_HAVE_DTAB :0)|
	    (pm->use_ptab   ?PFLAG_HAVE_PTAB :0)|
	    (pm->do_sel     ?PFLAG_DO_SEL    :0)|
	    (pm->fixed_len  ?PFLAG_DO_LEN    :0)|
	    (pm->do_dedup   ?PFLAG_DO_DEDUP  :0)|
	    (pm->store_qmap ?PFLAG_HAVE_QMAP :0);
    }

    for (i = gp->max_sel; i < 256; i++)
	gp->stab[i] = gp->stab[gp->max_sel-1];

    return 0;
}
#endif

#define NSYM 256
#define STEP 8
#include "htscodecs/c_simple_model.h"
#undef NSYM
#undef STEP

#define NSYM 4
#include "htscodecs/c_small_model.h"

#undef NSYM
#define NSYM 2
#include "htscodecs/c_small_model.h"

// An order-N arithmetic encoder, dedicated to sequence contexts.
char *encode_seq(unsigned char *in,  unsigned int in_size,
		 unsigned int *len, int nrecords, int both_strands,
		 int ctx_size, unsigned int *out_size) {
    char *out = malloc(in_size + 100);
    if (!out)
	return NULL;

    const int msize = 1<<(2*ctx_size);
    const int mask = msize-1;
    int i;

    SMALL_MODEL(4,_) *seq_model = malloc(msize * sizeof(*seq_model));

    for (i = 0; i < msize; i++)
	SMALL_MODEL(4,_init)(&seq_model[i]);

    SMALL_MODEL(2,_) state_model[3];
    SIMPLE_MODEL(256,_) run_len[3];
    for (i = 0; i < 3; i++) {
	SMALL_MODEL(2,_init)(&state_model[i]);
	SIMPLE_MODEL(256,_init)(&run_len[i], 256);
    }

    SIMPLE_MODEL(256,_) literal;
    SIMPLE_MODEL(256,_init)(&literal, 256);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out);
    RC_StartEncode(&rc);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    int last  = 0x007616c7 & mask;
    int last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask; // both strands mode

    int L[256];
    for (int i = 0; i < 256; i++)
	L[i] = 4; // N
    L['A'] = 0;
    L['C'] = 1;
    L['G'] = 2;
    L['T'] = 3;

    L['a'] = 0x80;
    L['c'] = 0x81;
    L['g'] = 0x82;
    L['t'] = 0x83;

    // Transition table to stored code:
    //    uc lc N
    // uc -  0  1
    // lc 0  -  1
    // N  0  1  -
    enum { uc_ACGT = 0, lc_ACGT = 1, other = 2 } state = uc_ACGT;

    int nseq = 0;
    int seq_len = len[nseq++];
    for (int i = 0; i < in_size; ) {//i++) {
	// Count size of consecutive symbols matching the same state
	int j, run, r2;
	switch (state) {
	case uc_ACGT: // uppercase ACGT symbols
	    for (j = i; j < in_size; j++) {
		if (L[in[j]] >= 4)
		    break;
	    }
	    break;

	case lc_ACGT: // lowercase acgt symbols
	    for (j = i; j < in_size; j++) {
		if (L[in[j]] < 0x80)
		    break;
	    }
	    break;

	case other: // ambiguity codes
	    for (j = i; j < in_size; j++) {
		if (L[in[j]] != 4)
		    break;
	    }
	    break;
	}

	// Encode the run length
	r2 = run = j-i;
	for (;;) {
	    //fprintf(stderr, "Encode %d of %d for state %d\n",
	    //    MIN(255, r2), run, state);
	    SIMPLE_MODEL(256, _encodeSymbol)(&run_len[state], &rc,
					     MIN(255, r2));
	    if (r2 >= 255)
		r2 -= 255;
	    else
		break;
	}

	// Encode the symbols
	switch (state) {
	case uc_ACGT:
	case lc_ACGT:
	    for (j = 0; j < run; j++) {
		unsigned char b = L[in[i+j]] & 3;
		SMALL_MODEL(4, _encodeSymbol)(&seq_model[last], &rc, b);

		last = ((last<<2) + b) & mask;
		int pf = ((last<<6)&mask)
			 +(i+j+3<in_size
			   ?L[in[i+j+1]]*16+L[in[i+j+2]]*4+L[in[i+j+3]]
			   :0);
#ifdef __SSE__
		_mm_prefetch((const char *)&seq_model[pf], _MM_HINT_T0);
#endif

		// 0.7% and 3.2% smaller for _.FQ and _.fq respectively
		// (at ctx_size 12), but 45% more CPU for seq encoding.
		if (both_strands) {
		    int b2 = last2 & 3;
		    last2 = last2/4 + ((3-b) << (2*ctx_size-2));
		    SMALL_MODEL(4, _updateSymbol)(&seq_model[last2], b2);

		    // ~25% speed gain by prefetching bottom strand too
		    uint32_t i3 = i+j+3 < in_size
			? L[in[i+j+1]] + L[in[i+j+2]]*4 + L[in[i+j+3]]*16
			: 0;
		    i3 = (0x3f - i3) << (2*ctx_size-6);
		    pf = i+j+3 < in_size ? (last2>>6) +i3 : 0;
		    _mm_prefetch((const char *)&seq_model[pf], _MM_HINT_T0);
		}

		// In theory we should reset context for each new sequence
		// as there is no obvious correlation between one sequence
		// and the next.  In practice the difference is 1-2%.
		// It only costs 1-2% CPU too, so worth doing.
		//
		//          -s5/100MB             -s7 -b / 512MB
		// Without: 76655812              60731188
		// With:    75789505 -1.1%        59638082 -1.8%
		// Slowdown 0.2%                  1.9%
		if (--seq_len == 0 && i+j+1 < in_size) {
		    if (nseq >= nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	    break;

	case other:
	    for (j = 0; j < run; j++) {
		SIMPLE_MODEL(256, _encodeSymbol)(&literal, &rc, in[i+j]);
		if (--seq_len == 0 && i+j+1 < in_size) {
		    if (nseq >= nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	}

	i += run;
	if (i >= in_size)
	    break;

	// Encode switch to next state
	switch(L[in[i]]) {
	case 0: case 1: case 2: case 3:
	    //fprintf(stderr, "state %d to 0, => 0\n", state);
	    SMALL_MODEL(2, _encodeSymbol)(&state_model[state], &rc, 0);
	    state = uc_ACGT;
	    break;
	case 0x80: case 0x81: case 0x82: case 0x83:
	    //fprintf(stderr, "state %d to 1, => %d\n", state, state==other);
	    SMALL_MODEL(2, _encodeSymbol)(&state_model[state], &rc,
					  state == other);
	    state = lc_ACGT;
	    break;
	default:
	    //fprintf(stderr, "state %d to 2, => 1\n", state);
	    SMALL_MODEL(2, _encodeSymbol)(&state_model[state], &rc, 1);
	    state = other;
	    break;
	}
    }

    RC_FinishEncode(&rc);
    *out_size = RC_OutSize(&rc);

    free(seq_model);

    return out;
}

char *decode_seq(unsigned char *in,  unsigned int in_size,
		 unsigned int *len, int nrecords, int both_strands,
		 int ctx_size, unsigned int out_size) {
    char *out = malloc(out_size);
    if (!out)
	return NULL;

    const int msize = 1<<(2*ctx_size);
    const int mask = msize-1;
    int i;

    SMALL_MODEL(4,_) *seq_model = malloc(msize * sizeof(*seq_model));

    // Do histogram to get observed values.
    // Then set m to max number of elements in histogram.
    for (i = 0; i < msize; i++)
	SMALL_MODEL(4,_init)(&seq_model[i]);

    SMALL_MODEL(2,_) state_model[3];
    SIMPLE_MODEL(256,_) run_len[3];
    for (i = 0; i < 3; i++) {
	SMALL_MODEL(2,_init)(&state_model[i]);
	SIMPLE_MODEL(256,_init)(&run_len[i], 256);
    }

    SIMPLE_MODEL(256,_) literal;
    SIMPLE_MODEL(256,_init)(&literal, 256);

    RangeCoder rc;
    RC_SetInput(&rc, (char *)in, (char *)in+in_size);
    RC_StartDecode(&rc);

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    int last  = 0x007616c7 & mask;
    int last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask; // both strands mode

    // Transition table to stored code:
    //    uc lc N
    // uc -  0  1
    // lc 0  -  1
    // N  0  1  -
    enum { uc_ACGT = 0, lc_ACGT = 1, other = 2 } state = uc_ACGT;

    int nseq = 0;
    int seq_len = len[nseq++];
    for (int i = 0; i < out_size;) {
	int j, run = 0, r2;
	// Fetch run length
	do {
	    r2 = SIMPLE_MODEL(256, _decodeSymbol)(&run_len[state], &rc);
	    run += r2;
	    //fprintf(stderr, "Decode %d of %d for state %d\n", r2, run, state);
	} while (r2 == 255);

	if (i + run > out_size)
	    // or error as it's malformed data
	    run = out_size - i;

	// Decode
	switch (state) {
	case uc_ACGT:
	case lc_ACGT: {
	    char *bases = state==lc_ACGT ? "acgt" : "ACGT";
	    for (j = 0; j < run; j++) {
		unsigned char b =
		    SMALL_MODEL(4, _decodeSymbol)(&seq_model[last], &rc);
		last = ((last<<2) + b) & mask;
#ifdef __SSE__
		_mm_prefetch((const char *)&seq_model[(last<<4)&mask],
			     _MM_HINT_T0);
#endif
		out[i+j] = bases[b];

		if (both_strands) {
		    int b2 = last2 & 3;
		    last2 = last2/4 + ((3-b) << (2*ctx_size-2));
		    SMALL_MODEL(4, _updateSymbol)(&seq_model[last2], b2);
	        }

		if (--seq_len == 0 && i+j+1 < out_size) {
		    if (nseq >=  nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	    break;
	}

	case other: // ambiguity codes
	    for (j = 0; j < run; j++) {
		out[i+j] = SIMPLE_MODEL(256, _decodeSymbol)(&literal, &rc);
		if (--seq_len == 0 && i+j+1 < out_size) {
		    if (nseq >=  nrecords) {
			free(out);
			free(seq_model);
			return NULL;
		    }
		    seq_len = len[nseq++];
		    last = 0x007616c7 & mask;
		    last2 = (0x2c6b62ff >> (32 - 2*ctx_size)) & mask;
		}
	    }
	    break;
	}

	i += run;
	if (i >= out_size)
	    break;

	// Next state
	int nstate = SMALL_MODEL(2, _decodeSymbol)(&state_model[state], &rc);
	switch (state) {
	case uc_ACGT:
	    state = nstate ? other : lc_ACGT;
	    break;
	case lc_ACGT:
	    state = nstate ? other : uc_ACGT;
	    break;
	case other:
	    state = nstate ? lc_ACGT : uc_ACGT;
	    break;
	}
    }

    RC_FinishDecode(&rc);

    free(seq_model);

    return out;
}

static char *encode_names(unsigned char *name_buf,  unsigned int name_len,
			  int strat, int level, unsigned int *out_size,
			  unsigned int *fq_flags, int num_records, int is_paired) {
    // TODO: work out a better maximum bound
    char *nout = malloc(name_len*2+1000), *cp = nout;
    if (!nout)
	return NULL;
    
    *(uint32_t *)cp = name_len; cp += 4;
    *cp++ = strat; // name method;

    if (strat == 0) {
	unsigned char *lzp_out = (unsigned char *)cp;
	unsigned int clen = lzp(name_buf, name_len, lzp_out);
	unsigned char *out = rans_compress_4x16(lzp_out, clen, &clen, 5);
	if (!out) {
	    free(nout);
	    return NULL;
	}
	*(uint32_t *)cp = clen; cp += 4;
	memcpy(cp, out, clen);  cp += clen;
	free(out);

    } else if (strat == 1) {
	int clen;
	unsigned char *out = tok3_encode_names((char *)name_buf, name_len,
					       level, 0, &clen, NULL);
	if (!out) {
	    free(nout);
	    return NULL;
	}
	*(uint32_t *)cp = clen; cp += 4;
	memcpy(cp, out, clen);  cp += clen;
	free(out);

    } else {
	char *n1 = malloc(name_len);
	if (!n1) {
	    free(nout);
	    return NULL;
	}
	char *n2 = malloc(name_len);
	if (!n2) {
	    free(n1);
	    free(nout);
	    return NULL;
	}
	unsigned char *flag = malloc(name_len/2);  // Worst case
	if (!flag) {
	    free(n1);
	    free(n2);
	    free(nout);
	    return NULL;
	}
	char *cp1 = n1, *cp2 = n2;
	int i = 0, nr = 0;
	// Flag bit 0: has "/NUM"
	// Flag bit 1: /1 vs /2
	// Flag bit 2: has a comment
	// Flag bit 3: space vs tab before comment
	while (i < name_len) {
	    int j, f = 0;
	    int w1end = 0;
	    int w2start = 0;
	    int w2end = 0;
	    for (j = i; j < name_len; j++) {
		if (name_buf[j] == '\0') {
		    w2end = j;
		    break;
		}
		if (!w2start && (name_buf[j] == ' ' ||
				 name_buf[j] == '\t')) {
		    w1end = j;
		    w2start = j+1;
		    f |= 4; // FLAG: has comment
		}
	    }

	    if (!w1end)
		w1end = j;
	    if (!w2end)
		w2end = j;

	    if (w2start)
		// FLAG: space vs tab
		f |= name_buf[w2start-1] == ' ' ? 0 : 8;

	    if (w1end>1 && name_buf[w1end-2] == '/') {
		// FLAG /1 or /2
		if (name_buf[w1end-1] == '1')
		    f |= 1, w1end -= 2;
		else if (name_buf[w1end-1] == '2')
		    f |= 3, w1end -= 2;
	    }

	    flag[nr++] = f;
	    memcpy(cp1, &name_buf[i], w1end-i);
	    cp1[w1end-i]=0;
	    cp1 += w1end-i+1;

	    if (w2start) {
		memcpy(cp2, &name_buf[w2start], w2end-w2start);
		cp2[w2end-w2start] = 0;
		cp2 += w2end-w2start+1;
	    }

	    i = j+1;
	}

	int clen1;
	unsigned int clenf, clen2 = 0;
	unsigned char *out = tok3_encode_names(n1, cp1-n1, level, 0, &clen1,
					       NULL);
	if (!out) {
	    free(n1);
	    free(n2);
	    free(flag);
	    free(nout);
	    return NULL;
	}
	unsigned char *outf = rans_compress_4x16(flag, nr, &clenf, 129);
	if (!outf) {
	    free(out);
	    free(n1);
	    free(n2);
	    free(flag);
	    free(nout);
	    return NULL;
	}
	unsigned char *out2 = NULL;
	if (cp2 != n2) {
	    unsigned char *lzp_out = malloc((cp2-n2)*2);
	    if (!lzp_out) {
		free(out);
		free(outf);
		free(n1);
		free(n2);
		free(flag);
		free(nout);
		return NULL;
	    }
	    clen2 = lzp((unsigned char *)n2, cp2-n2, lzp_out);
	    out2 = rans_compress_4x16(lzp_out, clen2, &clen2, 5);
	    free(lzp_out);
	    if (!out2) {
		free(out);
		free(outf);
		free(n1);
		free(n2);
		free(flag);
		free(nout);
		return NULL;
	    }
	}

	unsigned int clen = clen1 + clenf + clen2 + 8;

	*(uint32_t *)cp = clen;  cp += 4;
	*(uint32_t *)cp = clen1; cp += 4;
	*(uint32_t *)cp = clenf; cp += 4;

	memcpy(cp, out, clen1);  cp += clen1;
	free(out);

	memcpy(cp, outf, clenf); cp += clenf;
	free(outf);

	if (out2) {
	    memcpy(cp, out2, clen2); cp += clen2;
	    free(out2);
	}
	free(n1);
	free(n2);
	free(flag);
    }

    *out_size = cp - nout;
    return nout;
}

static char *decode_names(unsigned char *comp,  unsigned int c_len,
			  unsigned int u_len, int strat, unsigned int **out_flags, int *out_num_records) {
    unsigned char *out;

    if (strat == 0) {
	unsigned int ru_len;
	unsigned char *rout = rans_uncompress_4x16(comp, c_len, &ru_len);
	if (!rout)
	    goto err;
	out = malloc(u_len);
	if (!out) {
	    free(rout);
	    goto err;
	}
	u_len = unlzp(rout, ru_len, out);
	free(rout);
	// No flag array for strat 0
	if (out_flags)
	    *out_flags = NULL;
	if (out_num_records)
	    *out_num_records = 0;
    } else if (strat == 1) {
	out = tok3_decode_names(comp, c_len, &u_len);
	if (!out)
	    goto err;
	// No flag array for strat 1
	if (out_flags)
	    *out_flags = NULL;
	if (out_num_records)
	    *out_num_records = 0;
    } else {
	uint32_t clen1 = *(uint32_t *)comp;
	uint32_t clenf = *(uint32_t *)(comp+4);
	
	// Sanity check: ensure lengths are valid
	if (c_len < clen1 + clenf + 8) {
	    fprintf(stderr, "ERROR: Invalid compressed name data (c_len=%u, clen1=%u, clenf=%u)\n",
		    c_len, clen1, clenf);
	    goto err;
	}
	
	uint32_t clen2 = c_len - clen1 - clenf - 8;

	// Uncompress 3 separate components
	unsigned int u_len1, u_lenf, u_len2 = 0;
	unsigned char *out1 = tok3_decode_names(comp+8, clen1, &u_len1);
	if (!out1) {
	    fprintf(stderr, "ERROR: tok3_decode_names failed (clen1=%u, expected u_len=%u)\n",
		    clen1, u_len);
	    goto err;
	}
	unsigned char *outf = rans_uncompress_4x16(comp+8+clen1, clenf,
						   &u_lenf);
	if (!outf) {
	    fprintf(stderr, "ERROR: rans_uncompress_4x16 failed for flags (clenf=%u)\n",
		    clenf);
	    free(out1);
	    goto err;
	}
	unsigned char *out2 = NULL;
	if (clen2) {
	    unsigned int rulen;
	    unsigned char *rout = rans_uncompress_4x16(comp+8+clen1+clenf,
						       clen2, &rulen);
	    if (!rout) {
		fprintf(stderr, "ERROR: rans_uncompress_4x16 failed for comments (clen2=%u)\n",
			clen2);
		free(out1);
		free(outf);
		goto err;
	    }
	    out2 = malloc(u_len);
	    if (!out2) {
		fprintf(stderr, "ERROR: malloc failed for out2 (u_len=%u)\n", u_len);
		free(out1);
		free(outf);
		free(rout);
		goto err;
	    }
	    u_len2 = unlzp(rout, rulen, out2);
	    free(rout);
	}

	// Allocate array to store decoded flags
	unsigned int *decoded_flags = NULL;
	if (out_flags) {
	    decoded_flags = malloc(u_lenf * sizeof(unsigned int));
	    if (!decoded_flags) {
		free(out1);
		free(outf);
		free(out2);
		goto err;
	    }
	}

	// Stitch together ID + flag + comment
	unsigned char *cp1 = out1, *cp1_end = out1+u_len1;
	unsigned char *cpf = outf, *cpf_end = outf+u_lenf;
	unsigned char *cp2 = out2, *cp2_end = out2 + u_len2;
	
	// Allocate output buffer with extra space for potential /1 or /2 suffixes
	// Each record could add "/1" or "/2" (2 bytes), so add 2*u_lenf extra bytes
	// Check for potential overflow before calculations
	if (u_lenf > SIZE_MAX / 2) {
	    fprintf(stderr, "ERROR: Too many records for suffix calculation (u_lenf=%u)\n", u_lenf);
	    free(out1);
	    free(outf);
	    free(out2);
	    if (decoded_flags)
		free(decoded_flags);
	    goto err;
	}
	size_t suffix_space = (size_t)u_lenf * 2;
	size_t out_size = (size_t)u_len + suffix_space;
	// Sanity check: detect wraparound overflow
	// out_size < u_len means addition wrapped around
	if (out_size < u_len) {
	    fprintf(stderr, "ERROR: Output size overflow in decode_names (u_len=%u, u_lenf=%u)\n",
		    u_len, u_lenf);
	    free(out1);
	    free(outf);
	    free(out2);
	    if (decoded_flags)
		free(decoded_flags);
	    goto err;
	}
	out = malloc(out_size);
	if (!out) {
	    free(out1);
	    free(outf);
	    free(out2);
	    if (decoded_flags)
		free(decoded_flags);
	    goto err;
	}
	unsigned char *cp  = out,  *cp_end = out + out_size;  // Use out_size not u_len
	unsigned char *last_cp = NULL;
	int record_idx = 0;
	while (cp < cp_end) {
	    while (cp1 < cp1_end && cp < cp_end && *cp1)
		*cp++ = *cp1++;
	    cp1++;

	    int flag = 0;
	    if (cpf < cpf_end)
		flag = *cpf++;
	    // Need space for 2 bytes: '/' and '1'/'2'
	    if ((flag & 1) && cp + 1 < cp_end) {
		*cp++ = '/';
		*cp++ = (flag & 2) ? '2' : '1';
	    }
		
	    if (flag & 4 && cp < cp_end)
		*cp++ = (flag & 8) ? '\t' : ' ';

	    if (cp2) {
		while (cp2 < cp2_end && cp < cp_end && *cp2)
		    *cp++ = *cp2++;
		cp2++;
	    }

	    // Store the flag converted to FQZ_FREAD2 format
	    if (decoded_flags && record_idx < u_lenf) {
		// Convert name encoding flag to FQZ_FREAD2 flag
		// flag bit 0 = has /NUM, flag bit 1 = /2 (if set, otherwise /1)
		// For /2 reads: flag = 3 (bits 0 and 1 both set)
		// For /1 reads: flag = 1 (only bit 0 set)
		decoded_flags[record_idx] = ((flag & 3) == 3) ? FQZ_FREAD2 : 0;
	    }
	    record_idx++;

	    if (cp == last_cp)
		// ran out of data early; avoids looping forever
		break;

	    if (cp < cp_end) {
		*cp++ = 0;
	    } else {
		free(out1);
		free(outf);
		free(out2);
		if (decoded_flags) {
		    free(decoded_flags);
		    if (out_flags)
			*out_flags = NULL;
		}
		goto err;
	    }
	    last_cp = cp;
	}

	free(out1);
	free(outf);
	free(out2);
	
	if (out_flags)
	    *out_flags = decoded_flags;
	if (out_num_records)
	    *out_num_records = record_idx;
    }

    return (char *)out;

 err:
    free(out);
    return NULL;
}

// nstrat 0 = LZP+rans
// nstrat 1 = TOK3
// nstrat 2 = TOK3 name, LZP+rans comment
typedef struct {
    int nstrat, sstrat, qstrat;
    int scustom;                 // explicit -S -B user options
    int nlevel, slevel, qlevel;  // explicit level specified
    int nauto,  sauto,  qauto;   // generic -1 to -9 options (bit-wise levels)
    int verbose;
    int both_strands;
    uint32_t blk_size;
    int nthread;
    int plus_name;               // output name on third line (+name)
    int check_only;              // only verify integrity, don't decompress
    int inspect_only;            // inspect file and show metadata
    int verify_crc;              // verify CRC during decompression (enabled by default for v1.1)
    int paired_mode;             // 1 if processing paired-end reads (interleaved), 0 otherwise
} opts;

typedef struct {
    int64_t nblock;
    int64_t nusize, ncsize, ntime;
    int64_t susize, scsize, stime;
    int64_t qusize, qcsize, qtime;
    int64_t lusize, lcsize, ltime;
    int nmeth, smeth, qmeth, lmeth;
} timings;

static inline uint64_t tvdiff(struct timeval *tv1, struct timeval *tv2) {
    return (tv2->tv_sec - tv1->tv_sec) * 1000000
	+ tv2->tv_usec - tv1->tv_usec;
}

static inline void update_stats(timings *t,
		  int column, // 0=name 1=seq 2=qual 3=length
		  int64_t usize, int csize, int time) {
    switch (column) {
    case 0:
	t->nusize += usize;
	t->ncsize += csize;
	t->ntime  += time;
	break;
    case 1:
	t->susize += usize;
	t->scsize += csize;
	t->stime  += time;
	break;
    case 2:
	t->qusize += usize;
	t->qcsize += csize;
	t->qtime  += time;
	break;
    case 3:
	t->lusize += usize;
	t->lcsize += csize;
	t->ltime  += time;
	break;
    }
}

void append_timings(timings *t1, timings *t2, int verbose) {
    t1->nblock++;

    t1->nusize += t2->nusize;
    t1->ncsize += t2->ncsize;
    t1->ntime  += t2->ntime;

    t1->susize += t2->susize;
    t1->scsize += t2->scsize;
    t1->stime  += t2->stime;

    t1->qusize += t2->qusize;
    t1->qcsize += t2->qcsize;
    t1->qtime  += t2->qtime;

    t1->lusize += t2->lusize;
    t1->lcsize += t2->lcsize;
    t1->ltime  += t2->ltime;

    if (verbose) {
	fprintf(stderr, "Names   %11ld to %11ld in %.2f sec method %d\n",
		t2->nusize, t2->ncsize, t2->ntime/1e6, t2->nmeth);
	fprintf(stderr, "Lengths %11ld to %11ld in %.2f sec method %d\n",
		t2->lusize, t2->lcsize, t2->ltime/1e6, t2->lmeth);
	fprintf(stderr, "Seqs    %11ld to %11ld in %.2f sec method %d\n",
		t2->susize, t2->scsize, t2->stime/1e6, t2->smeth);
	fprintf(stderr, "Quals   %11ld to %11ld in %.2f sec method %d\n\n",
		t2->qusize, t2->qcsize, t2->qtime/1e6, t2->qmeth);
    }
}

#define APPEND_OUT(dat, len) do {			\
    if (*out_size < comp_sz + (len)) {			\
        *out_size = (comp_sz + (len))*1.5 + 1000;	\
	comp = realloc(comp, *out_size);		\
    }							\
    memcpy(comp+comp_sz, (dat), (len));			\
    comp_sz += (len);					\
} while(0);


// Updates the metrics counters and returns the methods to use.
// Method returned is a bitfield of (1<<method_num).
int metrics_method(int sec) {
    pthread_mutex_lock(&metric_m);

    if (stats[sec].review <= 0) {
	stats[sec].review = METRICS_REVIEW;
	stats[sec].trial  = METRICS_TRIAL;
	memset(stats[sec].usize, 0, M_LAST * sizeof(*stats[sec].usize));
	memset(stats[sec].csize, 0, M_LAST * sizeof(*stats[sec].csize));
	memset(stats[sec].count, 0, M_LAST * sizeof(*stats[sec].count));
    }

    int method;
    if (stats[sec].trial>0) {
	// Under evaluation => all methods used
	method = method_avail[sec];
    } else if (stats[sec].trial <= 0 && stats[sec].trial > -99999) {
	// Trial finished => select best method and set review timer

	int m, best_m = 0;
	//uint32_t best_sz = UINT_MAX;
	double best_sz = 1e30;
	for (m = 0; m < M_LAST; m++) {
	    // TODO: parameterise by speed too, plus small block offset?
	    if (stats[sec].usize[m] && best_sz > (stats[sec].csize[m]+1.0)/stats[sec].usize[m]) {
		best_sz = (stats[sec].csize[m]+1.0)/stats[sec].usize[m];
		best_m = m;
	    }
	}
	//fprintf(stderr, "Choose best method %d for sec %d\n", best_m, sec);
	method_used[sec] = best_m;
	method = 1<<method_used[sec];

	stats[sec].trial  = -99999;

	// TODO: methods that consistently get rejected can be removed from
	// the method_avail.  This is a second level based on
	// total accumulated size stats.
    } else {
	// Repeat best method until review counter hits zero again
	stats[sec].review--;
	method = 1<<method_used[sec];
    }

    pthread_mutex_unlock(&metric_m);

    return method;
}

// Update the metrics for a given section and method.
// Method parameter isn't a bitfield, but the method number itself.
void metrics_update(int sec, int method, int64_t usize, int64_t csize) {
    if (stats[sec].trial <= 0)
	return; // done
    //pthread_mutex_lock(&metric_m);
    stats[sec].usize[method] += usize;
    stats[sec].csize[method] += csize;
    stats[sec].count[method]++;
    //pthread_mutex_unlock(&metric_m);
    //fprintf(stderr, "Section %d  method %d  size %ld\n", sec, method, csize);
}

// TODO: return buffer + meta, so we don't need memcpy and memmove calls.
char *compress_with_methods(fqz_gparams *gp,  opts *arg, fastq *fq,
			    uint32_t methods,
			    int sec, char *in, unsigned int in_size,
			    unsigned int *out_size, int *strat,
			    int *meth_used) {
    uint8_t *best_comp = NULL;
    uint32_t best_sz = UINT_MAX;
    int      best_strat = 0, best_method = 0;
    char *out;
    int m;
    size_t out_len;

    pthread_mutex_lock(&metric_m);
    int in_trial = stats[sec].trial > 0;
    pthread_mutex_unlock(&metric_m);

    metrics local_stats = {{0}};

    for (m = 0; m < M_LAST; m++) {
	if (!(methods & (1<<m)))
	    continue;

	out_len = UINT_MAX;

	switch (m) {
	case RANS0:
	case RANS1:
	case RANS64:
	case RANS65:
	case RANS128:
	case RANS129:
	case RANS192:
	case RANS193: {
	    int order[] = {0,1,64,65,128,129,192,193};
	    *strat = 0;
	    out = (char *)rans_compress_4x16((uint8_t *)in, in_size,
					     out_size, order[m-RANS0]);
	    out_len = *out_size;
	    break;
	}

	case RANSXN1:
	    if (!fq->fixed_len) {
		out = NULL;
		break;
	    }
	    *strat = 0;
	    out = (char *)rans_compress_4x16((uint8_t *)in, in_size, out_size,
					     (fq->fixed_len<<8)+9);
	    out_len = *out_size;
	    break;

	case LZP3: {
	    unsigned char *lzp_out = malloc(in_size*2+1000);
	    unsigned int lzp_len = lzp(in, in_size, lzp_out);
	    out = (char *)rans_compress_4x16(lzp_out, lzp_len, out_size, 5);
	    free(lzp_out);
	    out_len = *out_size;
	    *strat = LZP3;
	    break;
	}

	case TLZP3:
	    out = encode_names(in, in_size, 0 /* LZP + rANS o5 */,
			       (m-TOK3_3)*2+3, out_size, fq->flag, fq->num_records, arg->paired_mode);
	    out_len = *out_size;
	    break;

	case TOK3_3:
	case TOK3_5:
	case TOK3_7:
	case TOK3_9:
	    out = encode_names(in, in_size, 1 /* TOK3 */,
			       (m-TOK3_3)*2+3, out_size, fq->flag, fq->num_records, arg->paired_mode);
	    out_len = *out_size;
	    break;

	case TOK3_3_LZP:
	case TOK3_5_LZP:
	case TOK3_7_LZP:
	case TOK3_9_LZP:
	    out = encode_names(in, in_size, 2 /* TOK3+LZP */,
			       (m-TOK3_3_LZP)*2+3, out_size, fq->flag, fq->num_records, arg->paired_mode);
	    out_len = *out_size;
	    break;

	case SEQ10:
	case SEQ12:
	case SEQ12B:
	case SEQ13B:
	case SEQ14B: {
	    int slevel[]   = {10,12,12,13,14};
	    int both_str[] = {0, 0, 1, 1, 1};

	    int s = m-SEQ10;

	    *strat = (slevel[s]<<4) | (both_str[s]<<3) | 1;
	    out = encode_seq(in, in_size, fq->len, fq->num_records,
			     both_str[s], slevel[s], out_size);
	    out_len = *out_size;
	    break;
	}

	case SEQ_CUSTOM:
	    *strat = (arg->slevel<<4) | (arg->both_strands<<3) | 1;
	    out = encode_seq(in, in_size, fq->len, fq->num_records,
			     arg->both_strands, arg->slevel, out_size);
	    out_len = *out_size;
	    break;

	case FQZ0:
	case FQZ1:
	case FQZ2:
	case FQZ3:
	case FQZ4: {
	    *strat = 1;
	    fqz_slice *s = malloc(fq->num_records * sizeof(*s));
	    s->num_records = fq->num_records;
	    s->len = fq->len;
	    s->flags = fq->flag;
	    s->seq = malloc(fq->num_records * sizeof(char *));
	    int i, j;
	    for (i = j = 0; i < fq->num_records; j += fq->len[i++])
		s->seq[i] = (unsigned char *)fq->seq_buf + j;

	    // FIXME: expose fqz_pick_parameters function so we
	    // can initialise it here and then also turn off DO_LEN.
	    out = fqz_compress(4, s, in, in_size, 
			       &out_len, m-FQZ0, gp);
	    free(s->seq);
	    free(s);
	    break;
	}

	default:
	    fprintf(stderr, "Unsupported method %d (set 0x%x) for section %d\n", m, methods, sec);
	    abort();
	}

	if (arg->verbose > 2) {
	    char *secstr[] = {"name", "length", "sequence", "quality"};
	    fprintf(stderr, "Try      %8s with method %2d %10d to %10d "
		    "bytes\n", secstr[sec], m, in_size, (uint32_t)out_len);
	}

	if (best_sz > out_len) {
	    best_sz = out_len;
	    best_method = m;
	    if (best_comp)
		free(best_comp);
	    best_comp = out;
	    best_strat = *strat;
	} else {
	    free(out);
	}

	local_stats.usize[m] = in_size;
	local_stats.csize[m] = out_len;
    }

    if (in_trial) {
	pthread_mutex_lock(&metric_m);
	for (m = 0; m < M_LAST; m++) {
	    if (!(methods & (1<<m)))
		continue;
	    metrics_update(sec, m, local_stats.usize[m], local_stats.csize[m]);
	}
	stats[sec].trial--;
	pthread_mutex_unlock(&metric_m);
    }

    out = best_comp;

    if (arg->verbose > 1) {
	char *secstr[] = {"name", "length", "sequence", "quality"};
	fprintf(stderr, "Compress %8s with method %2d %10d to %10d "
		"bytes\n", secstr[sec], best_method, in_size, best_sz);
    }

    *out_size = best_sz;
    *strat    = best_strat;
    *meth_used = best_method;
    return out;
}

// Encodes a single block of data
char *encode_block(fqz_gparams *gp, opts *arg, fastq *fq, timings *t,
		   unsigned int *out_size) {
    struct timeval tv1, tv2;

    // Starting guess
    *out_size = 1000;//arg->blk_size/4 + 10000;
    char *comp = malloc(*out_size), *out;
    unsigned int clen, comp_sz = 0;
    int strat = 0, method;

    // Reserve space for block size (will be filled at the end)
    uint32_t block_size_offset = comp_sz;
    comp_sz += 4;  // Reserve 4 bytes for block size
    
    APPEND_OUT(&fq->num_records, 4);
    
    // Reserve space for CRC32 (will be filled at the end)
    uint32_t crc_offset = comp_sz;
    comp_sz += 4;  // Reserve 4 bytes for CRC

    //----------
    // Names: tok3
    // Strat 0 = LZP + rANS
    // Strat 1 = Tok3
    // Strat 2 = Name(tok3)+Flag(RC)+Comment(LZP+rANS)
    gettimeofday(&tv1, NULL);
#if 1
    method = metrics_method(SEC_NAME);
    out = compress_with_methods(gp, arg, fq, method, SEC_NAME,
				fq->name_buf, fq->name_len,
				&clen, &strat, &t->nmeth);
#else
    out = encode_names((uint8_t *)fq->name_buf, fq->name_len, arg->nstrat,
		       arg->nlevel, &clen, fq->flag, fq->num_records, arg->paired_mode);
#endif
    APPEND_OUT(out, clen);
    free(out);

    gettimeofday(&tv2, NULL);
    update_stats(t, 0, fq->name_len, clen, tvdiff(&tv1, &tv2));

    //----------
    // Read lengths
    if (fq->fixed_len) {
	// Fixed length, with next byte holding the size of length
	unsigned char buf[5], nb = 1;
	nb += var_put_u32(buf+1, NULL, fq->fixed_len);
	buf[0] = nb-1;
	update_stats(t, 3, 4*fq->num_records, nb, 0);
	APPEND_OUT(buf, nb);
	t->lmeth = 1;
    } else {
	// Variable length (next byte 0), with 4 byte len followed
	// by var-int lengths.
	int i, nb = 0;
	unsigned char *buf = malloc(fq->num_records*5+5);
	buf[nb++] = 0;
	nb += 4; // comp.size placeholder

	for (i = 0; i < fq->num_records; i++)
	    nb += var_put_u32(buf+nb, NULL, fq->len[i]);
	*(uint32_t *)(buf+1) = nb-5;
	APPEND_OUT(buf, nb);
	free(buf);

	update_stats(t, 3, 4*fq->num_records, nb, 0);
	t->lmeth = 0;
    }

    //----------
    // Seq: rans or statistical modelling
    gettimeofday(&tv1, NULL);
    uint8_t  meta[9];

    method = metrics_method(SEC_SEQ);
    out = compress_with_methods(gp, arg, fq, method, SEC_SEQ,
				fq->seq_buf, fq->seq_len,
				&clen, &strat, &t->smeth);
    meta[0] = strat;

    *(uint32_t *)(&meta[1]) = fq->seq_len;
    *(uint32_t *)(&meta[5]) = clen;
    APPEND_OUT(meta, 9);
    APPEND_OUT(out, clen);
    free(out);

    gettimeofday(&tv2, NULL);
    update_stats(t, 1, fq->seq_len, clen+9, tvdiff(&tv1, &tv2));

    //----------
    // Qual: rans or fqz (skip for FASTA)
    if (!fq->is_fasta) {
        // FASTQ format - compress quality
        gettimeofday(&tv1, NULL);
        size_t out_len = 0;

        method = metrics_method(SEC_QUAL);
        out = compress_with_methods(gp, arg, fq, method, SEC_QUAL,
				    fq->qual_buf, fq->qual_len,
				    &clen, &strat, &t->qmeth);
        meta[0] = strat;
        //fprintf(stderr, "Qual %d -> %d via %d\n", fq->qual_len, clen, strat);

        *(uint32_t *)(&meta[1]) = fq->qual_len;
        *(uint32_t *)(&meta[5]) = clen;
        APPEND_OUT(meta, 9);
        APPEND_OUT(out, clen);
        free(out);

        gettimeofday(&tv2, NULL);
        update_stats(t, 2, fq->qual_len, clen+9, tvdiff(&tv1, &tv2));
    } else {
        // FASTA format - no quality
        meta[0] = 0;
        *(uint32_t *)(&meta[1]) = 0;
        *(uint32_t *)(&meta[5]) = 0;
        APPEND_OUT(meta, 9);
    }

    // Compute CRC32 for the block data (from position 12 onwards)
    // CRC is computed on all data after the CRC field itself
    uint32_t block_crc = crc32(0L, Z_NULL, 0);
    block_crc = crc32(block_crc, (unsigned char *)(comp + 12), comp_sz - 12);
    
    // Insert CRC at the reserved position (position 8, after block_size and num_records)
    *(uint32_t *)(comp + crc_offset) = block_crc;
    
    // Update block size to include CRC and all the data
    *(uint32_t *)(comp + block_size_offset) = comp_sz - 4;

    *out_size = comp_sz;

    return comp;
}

#define GET(ptr, len)				\
    do {					\
        if (in_off + (len) > in_size)		\
	    goto err;				\
	memcpy((ptr), in + in_off, (len));	\
	in_off += (len);			\
    } while(0);

fastq *decode_block(unsigned char *in, unsigned int in_size, timings *t, int file_version) {
    unsigned char *in_end = in + in_size, *comp, *out;
    uint32_t in_off = 0, nr;
    int i, j, err = 0;
    uint32_t u_len, c_len;
    uint8_t c;
    struct timeval tv1, tv2;
    uint32_t block_size;
    uint32_t block_crc_stored = 0, block_crc_computed = 0;
    
    // Read block size (new format)
    GET(&block_size, 4);
    
    GET(&nr, 4);
    
    // Read and verify CRC if this is v1.1 format
    if (file_version == 0) {  // v1.1 with CRC
        GET(&block_crc_stored, 4);
        
        // Compute CRC of the block data (from current position onwards)
        block_crc_computed = crc32(0L, Z_NULL, 0);
        block_crc_computed = crc32(block_crc_computed, in + in_off, block_size - 8);  // -8 for nr and crc fields
        
        if (block_crc_stored != block_crc_computed) {
            fprintf(stderr, "ERROR: Block CRC mismatch! File may be corrupted.\n");
            fprintf(stderr, "  Expected: 0x%08x, Got: 0x%08x\n", block_crc_stored, block_crc_computed);
            return NULL;
        }
    }
    
    fastq *fq = fastq_alloc(nr);

    // ----------
    // Name
    gettimeofday(&tv1, NULL);
    GET(&u_len, 4);
    GET(&c, 1);     // strategy
    GET(&c_len, 4);

    comp = in+in_off;
    in_off += c_len;

    unsigned int *decoded_flags = NULL;
    int num_decoded_records = 0;
    out = (unsigned char *)decode_names(comp, c_len, u_len, c, &decoded_flags, &num_decoded_records);
    if (!out) {
        fprintf(stderr, "ERROR: Failed to decode names\n");
        goto err;
    }
    fq->name_buf = (char *)out;
    fq->name_len = u_len;

    // Populate name indices and flags from decompressed names
    {
        char *np = fq->name_buf;
        int name_i = 0;
        int last_name = -1;
        for (i = 0; i < nr; i++) {
            fq->name[i] = name_i;
            
            // Use decoded flags if available (from strat==2), otherwise parse from name
            int flag = 0;
            if (decoded_flags && i < num_decoded_records) {
                flag = decoded_flags[i];
            } else {
                // Fallback: determine flag by checking name suffix or comparing with previous name
                int name_len = strlen(np);
                if (name_len > 1 && np[name_len-1] == '2' && np[name_len-2] == '/')
                    flag = FQZ_FREAD2;
                else if (last_name >= 0 &&
                    strcmp(fq->name_buf + fq->name[i], fq->name_buf + last_name) == 0)
                    flag = FQZ_FREAD2;
            }
            fq->flag[i] = flag;
            
            if (!flag)
                last_name = fq->name[i];
            
            name_i += strlen(np) + 1;  // +1 for null terminator
            np += strlen(np) + 1;
        }
        
        // Free decoded flags array
        free(decoded_flags);
    }

    gettimeofday(&tv2, NULL);
    t->ncsize += u_len;
    t->nusize += c_len;
    t->ntime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
    t->ntime += tv2.tv_usec - tv1.tv_usec;

    // ----------
    // Lengths
    GET(&c, 1);     // strategy, but also length.  Needed as len?
    if (c > 0) {
	// Fixed length
	uint32_t len;
	int vl;
	in_off += (vl = var_get_u32(in+in_off, in_end, &len));
	err |= vl == 0;

	for (i = 0; i < nr; i++)
	    fq->len[i] = len;
	t->lcsize += nr*4;
	t->lusize += c;
    } else {
	// Variable length
	uint32_t blen;
	GET(&blen, 4); // needed?  Doesn't seem it now
	for (i = 0; i < nr; i++) {
	    int vl;
	    in_off += (vl = var_get_u32(in+in_off, in_end, &fq->len[i]));
	    err |= vl == 0;
	}

	t->lcsize += nr*4;
	t->lusize += blen+5;
    }
    if (err)
	goto err;

    // ----------
    // Seq
    gettimeofday(&tv1, NULL);
    GET(&c, 1);
    GET(&u_len, 4);
    GET(&c_len, 4);
    comp = in+in_off;
    in_off += c_len;

    int slevel = c>>4;
    int both_strands = (c >> 3) & 1;

    if ((c & 7) == 1) {
	out = (uint8_t *)decode_seq(comp, c_len, fq->len, fq->num_records,
				    both_strands, slevel, u_len);
	if (!out) {
	    fprintf(stderr, "ERROR: Failed to decode sequence data (FQZ)\n");
	    goto err;
	}
    } else if (c == LZP3) {
	unsigned int ru_len;
	unsigned char *rout = rans_uncompress_4x16(comp, c_len, &ru_len);
	if (!rout) {
	    fprintf(stderr, "ERROR: Failed to decompress sequence data (rANS)\n");
	    goto err;
	}
	out = malloc(u_len);
	if (!out) {
	    fprintf(stderr, "ERROR: Failed to allocate memory for sequence\n");
	    free(rout);
	    goto err;
	}
	u_len = unlzp(rout, ru_len, out);
	free(rout);
    } else if (c == 0) {
	out = rans_uncompress_4x16(comp, c_len, &u_len);
	if (!out) {
	    fprintf(stderr, "ERROR: Failed to decompress sequence data (rANS)\n");
	    goto err;
	}
    } else {
	fprintf(stderr, "Unrecognised sequence strategy %d\n", c);
	goto err;
    }

    fq->seq_buf = (char *)out;
    fq->seq_len = u_len;
    for (i = 0; i < nr; i++)
	fq->seq[i] = i ? fq->seq[i-1] + fq->len[i-1] : 0;

    gettimeofday(&tv2, NULL);
    t->stime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
    t->stime += tv2.tv_usec - tv1.tv_usec;
    t->scsize += u_len;
    t->susize += c_len;

    // ----------
    // Qual (handle FASTA)
    gettimeofday(&tv1, NULL);
    GET(&c, 1);
    GET(&u_len, 4);
    GET(&c_len, 4);
    
    size_t out_len = 0;
    
    if (u_len == 0 && c_len == 0) {
        // FASTA format
        fq->is_fasta = 1;
        fq->qual_buf = calloc(1, 1);
        fq->qual_len = 0;
        for (i = 0; i < nr; i++)
            fq->qual[i] = 0;
    } else {
        // FASTQ format
        comp = in+in_off;
        in_off += c_len;
	    if (c == 0) {
		// Rans
		out = rans_uncompress_4x16(comp, c_len, &u_len);
		if (!out) {
		    fprintf(stderr, "ERROR: Failed to decompress quality data (rANS)\n");
		    goto err;
		}
		fq->qual_buf = (char *)out;
		fq->qual_len = out_len = u_len;
	    } else {
		// FQZComp qual
		fqz_slice s;
		s.num_records = fq->num_records;
		s.len = fq->len;
		s.flags = fq->flag;
		//s.seq = (unsigned char **)fq->seq;
		s.seq = (unsigned char **)malloc(fq->num_records * sizeof(char *));
		if (!s.seq) {
		    fprintf(stderr, "ERROR: Failed to allocate memory for sequence pointers\n");
		    goto err;
		}
		for (i = j = 0; i < fq->num_records; j += fq->len[i++])
		    s.seq[i] = (unsigned char *)fq->seq_buf + j;

		// pass lengths as NULL and fix fqz_decompress to cope?
		int *lengths = malloc(nr * sizeof(lengths));
		if (!lengths) {
		    fprintf(stderr, "ERROR: Failed to allocate memory for lengths\n");
		    free(s.seq);
		    goto err;
		}
		out = (uint8_t *)fqz_decompress((char *)comp, c_len, &out_len,
						lengths, nr, &s);
		if (!out) {
		    fprintf(stderr, "ERROR: Failed to decompress quality data (FQZ)\n");
		    free(s.seq);
		    free(lengths);
		    goto err;
		}
		fq->qual_buf = (char *)out;
		fq->qual_len = out_len;
		free(s.seq);
		free(lengths);
	    }
	for (i = 0; i < fq->qual_len; i++)
	    fq->qual_buf[i] += 33;
    }

    gettimeofday(&tv2, NULL);
    t->qtime += (tv2.tv_sec - tv1.tv_sec) * 1000000;
    t->qtime += tv2.tv_usec - tv1.tv_usec;
    t->qcsize += out_len;
    t->qusize += c_len;

    return fq;

 err:
    free(fq);
    return NULL;
}

typedef struct {
    fqz_gparams *gp;
    opts *arg;
    timings t;
    fastq *fq;
    char *comp;
    uint32_t clen;
    uint32_t usize;      // Uncompressed size (for index)
    uint32_t nrecords;   // Number of records (for index)
    int eof;
    int file_version;    // 0=v1.1 (with CRC), 1=v1.0 (no CRC), 2=old format
} enc_dec_job;

// Write FQZ5 file header
static int write_header(FILE *fp) {
    // Write magic number and version
    if (fwrite(FQZ5_MAGIC, 1, FQZ5_MAGIC_LEN, fp) != FQZ5_MAGIC_LEN)
        return -1;
    
    // Write index offset (0 for now, will be updated later)
    uint64_t index_offset = 0;
    if (fwrite(&index_offset, 1, 8, fp) != 8)
        return -1;
    
    return 0;
}

// Read FQZ5 file header
// Returns: 0 = v1.1 format with CRC, 1 = v1.0 format with CRC, 2 = old format no header, -1 = error
static int read_header(FILE *fp, uint64_t *index_offset) {
    char magic[FQZ5_MAGIC_LEN];
    
    if (fread(magic, 1, FQZ5_MAGIC_LEN, fp) != FQZ5_MAGIC_LEN)
        return -1;
    
    // Check for v1.1 format (with CRC)
    if (memcmp(magic, FQZ5_MAGIC, FQZ5_MAGIC_LEN) == 0) {
        if (fread(index_offset, 1, 8, fp) != 8)
            return -1;
        return 0; // v1.1 format with CRC
    }
    
    // Check for v1.0 format (no CRC)
    if (memcmp(magic, FQZ5_MAGIC_V10, FQZ5_MAGIC_LEN) == 0) {
        if (fread(index_offset, 1, 8, fp) != 8)
            return -1;
        return 1; // v1.0 format without CRC
    }
    
    // Not a FQZ5 file - old format without header
    // For backward compatibility, rewind to beginning
    fseek(fp, 0, SEEK_SET);
    *index_offset = 0;
    return 2; // Old format
}

// Write index at end of file
static int write_index(FILE *fp, fqz5_index *idx) {
    if (!idx || idx->nblocks == 0)
        return 0; // Nothing to write
    
    // Write index magic
    if (fwrite(FQZ5_INDEX_MAGIC, 1, FQZ5_INDEX_MAGIC_LEN, fp) != FQZ5_INDEX_MAGIC_LEN)
        return -1;
    
    // Write number of blocks
    if (fwrite(&idx->nblocks, 1, 4, fp) != 4)
        return -1;
    
    // Write index entries
    for (uint32_t i = 0; i < idx->nblocks; i++) {
        if (fwrite(&idx->entries[i].offset, 1, 8, fp) != 8)
            return -1;
        if (fwrite(&idx->entries[i].usize, 1, 4, fp) != 4)
            return -1;
        if (fwrite(&idx->entries[i].nrecords, 1, 4, fp) != 4)
            return -1;
    }
    
    return 0;
}

// Read index from file
static fqz5_index *read_index(FILE *fp, uint64_t index_offset) {
    if (index_offset == 0)
        return NULL; // No index
    
    if (fseek(fp, index_offset, SEEK_SET) != 0)
        return NULL;
    
    char magic[FQZ5_INDEX_MAGIC_LEN];
    if (fread(magic, 1, FQZ5_INDEX_MAGIC_LEN, fp) != FQZ5_INDEX_MAGIC_LEN)
        return NULL;
    
    if (memcmp(magic, FQZ5_INDEX_MAGIC, FQZ5_INDEX_MAGIC_LEN) != 0)
        return NULL;
    
    fqz5_index *idx = calloc(1, sizeof(*idx));
    if (!idx)
        return NULL;
    
    if (fread(&idx->nblocks, 1, 4, fp) != 4) {
        free(idx);
        return NULL;
    }
    
    idx->entries = calloc(idx->nblocks, sizeof(*idx->entries));
    if (!idx->entries) {
        free(idx);
        return NULL;
    }
    
    for (uint32_t i = 0; i < idx->nblocks; i++) {
        if (fread(&idx->entries[i].offset, 1, 8, fp) != 8 ||
            fread(&idx->entries[i].usize, 1, 4, fp) != 4 ||
            fread(&idx->entries[i].nrecords, 1, 4, fp) != 4) {
            free(idx->entries);
            free(idx);
            return NULL;
        }
    }
    
    return idx;
}

// Free index
static void free_index(fqz5_index *idx) {
    if (!idx)
        return;
    free(idx->entries);
    free(idx);
}

// Write trailer with overall file CRC
static int write_trailer(FILE *fp, uint32_t overall_crc, uint32_t nblocks) {
    // Write trailer magic
    if (fwrite(FQZ5_TRAILER_MAGIC, 1, FQZ5_TRAILER_MAGIC_LEN, fp) != FQZ5_TRAILER_MAGIC_LEN)
        return -1;
    
    // Write overall CRC
    if (fwrite(&overall_crc, 1, 4, fp) != 4)
        return -1;
    
    // Write number of blocks for verification
    if (fwrite(&nblocks, 1, 4, fp) != 4)
        return -1;
    
    return 0;
}

// Read and verify trailer
static int read_trailer(FILE *fp, uint32_t *overall_crc, uint32_t *nblocks) {
    char magic[FQZ5_TRAILER_MAGIC_LEN];
    
    // Try to read trailer magic
    if (fread(magic, 1, FQZ5_TRAILER_MAGIC_LEN, fp) != FQZ5_TRAILER_MAGIC_LEN)
        return -1; // No trailer (probably v1.0 file)
    
    if (memcmp(magic, FQZ5_TRAILER_MAGIC, FQZ5_TRAILER_MAGIC_LEN) != 0)
        return -1; // Not a valid trailer
    
    // Read overall CRC
    if (fread(overall_crc, 1, 4, fp) != 4)
        return -1;
    
    // Read number of blocks
    if (fread(nblocks, 1, 4, fp) != 4)
        return -1;
    
    return 0;
}

static void *encode_thread(void *arg) {
    enc_dec_job *j = (enc_dec_job *)arg;
    if (j->eof)
	return j;

    j->comp = encode_block(j->gp, j->arg, j->fq, &j->t, &j->clen);
    // j->usize and j->nrecords already set before dispatch
    fastq_free(j->fq);
    return j;
}

#define THREADED

// TODO: use async read and write threads too so main doesn't block on I/O
int encode(FILE *in_fp, FILE *out_fp, fqz_gparams *gp, opts *arg,
	   timings *t) {
    int rans_methods = (1<<RANS0) | (1<<RANS1) | (1<<RANS129) | (1<<RANS193);

    // Write file header with magic number
    if (write_header(out_fp) < 0)
        return -1;
    
    // Initialize index
    fqz5_index idx = {0};
    uint32_t idx_capacity = 1000;
    idx.entries = calloc(idx_capacity, sizeof(*idx.entries));
    if (!idx.entries)
        return -1;

    // Name
    if (arg->nauto) {
	method_avail[SEC_NAME] = arg->nauto;
    } else {
	if (arg->nstrat == 1)
	    method_avail[SEC_NAME] |= 1<<(TOK3_3 + arg->nlevel/2-1);
	else if (arg->nstrat == 2)
	    method_avail[SEC_NAME] |= 1<<(TOK3_3_LZP + arg->nlevel/2-1);
	else
	    method_avail[SEC_NAME] = 1<<TLZP3;
    }
    
    // Seq
    if (arg->scustom) {
	method_avail[SEC_SEQ] = 1<<SEQ_CUSTOM;
    } else {
	if (arg->sauto)
	    method_avail[SEC_SEQ] = arg->sauto;
	else if (arg->sstrat == 1)
	    method_avail[SEC_SEQ] = 1<<SEQ_CUSTOM;

	if (!method_avail[SEC_SEQ])
	    method_avail[SEC_SEQ] = rans_methods;
    }

    // Qual
    if (arg->qauto) {
	method_avail[SEC_QUAL] = arg->qauto;
    } else {
	if (arg->qstrat == 1) {
	    if (arg->qlevel == 4)
		method_avail[SEC_QUAL] = FQZ4;
	    else if (arg->qlevel == 3)
		method_avail[SEC_QUAL] = FQZ3;
	    else if (arg->qlevel == 1)
		method_avail[SEC_QUAL] = FQZ1;
	    else if (arg->qlevel == 2)
		method_avail[SEC_QUAL] = FQZ2;
	    else
		method_avail[SEC_QUAL] = FQZ0;
	} else {
	    method_avail[SEC_QUAL] = rans_methods;
	}
    }

#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n, 0);
    hts_tpool_result *r;
    enc_dec_job *j, *jr;
#endif

    char *in = malloc(arg->blk_size);
    int in_rem = 0;

    for(;;) {
	int nbytes = fread(in+in_rem, 1, arg->blk_size - in_rem, in_fp);
	if (nbytes < 0)
	    return -1;
	if (nbytes == 0)
	    break;
	nbytes += in_rem;

	fastq *fq = load_seqs(in, nbytes, &in_rem);
	memmove(in, in+in_rem, nbytes - in_rem);
	in_rem = nbytes - in_rem;

	if (!fq)
	    return -1;

	if (!fq->num_records) {
	    fastq_free(fq);
	    break;
	}

#ifdef THREADED
	// Dispatch a job
	j = calloc(1, sizeof(*j));
	j->gp = gp;
	j->arg = arg;
	memset(&j->t, 0, sizeof(j->t));
	j->fq = fq;
	j->eof = 0;
	j->usize = fq->seq_len;  // Total uncompressed bases
	j->nrecords = fq->num_records;

	int ret = -1;
	while (ret == -1) {
	    // Always dispatch, going over-size on queue
	    if ((ret=hts_tpool_dispatch2(p, q, encode_thread, j, -1)) != 0) {
		if (errno != EAGAIN)
		    goto err;
		//fprintf(stderr, "Couldn't dispatch\n");
	    }
	    //fprintf(stderr, "r=%d\n", ret);

	    // Check for a result.
	    // If input queue is oversize then do this blocking so we
	    // don't grow input queue indefinitely.
	    // TODO: same issue exists on decoder side.
	    do {
		if (hts_tpool_dispatch_would_block(p, q)) {
		    //fprintf(stderr, "would block: in=%d\n", q->n_input);
		    r = hts_tpool_next_result_wait(q);
		} else {
		    //fprintf(stderr, "would pass : in=%d\n", q->n_input);
		    r = hts_tpool_next_result(q);
		}
		if (r) {
		    jr = hts_tpool_result_data(r);
		    if (jr->eof) {
			end = 1;
		    } else {
			append_timings(t, &jr->t, arg->verbose);
			
			// Track block offset for index
			if (idx.nblocks >= idx_capacity) {
			    idx_capacity *= 2;
			    index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
			    if (!new_entries) {
				free(idx.entries);
				goto err;
			    }
			    idx.entries = new_entries;
			}
			idx.entries[idx.nblocks].offset = ftell(out_fp);
			idx.entries[idx.nblocks].usize = jr->usize;
			idx.entries[idx.nblocks].nrecords = jr->nrecords;
			idx.nblocks++;
			
			// Block now includes block_size at the start, no need for extra length
			fwrite(jr->comp, 1, jr->clen, out_fp);
			free(jr->comp);
		    }
		    hts_tpool_delete_result(r, 1);
		}
	    } while (r && hts_tpool_dispatch_would_block(p, q)); // necessary?
	}
#else
	uint32_t clen;
	t->nblock++;
	
	// Track block offset for index
	if (idx.nblocks >= idx_capacity) {
	    idx_capacity *= 2;
	    index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
	    if (!new_entries) {
		free(idx.entries);
		return -1;
	    }
	    idx.entries = new_entries;
	}
	idx.entries[idx.nblocks].offset = ftell(out_fp);
	idx.entries[idx.nblocks].usize = fq->seq_len;
	idx.entries[idx.nblocks].nrecords = fq->num_records;
	idx.nblocks++;
	
	char *out = encode_block(gp, arg, fq, t, &clen);
	// Block now includes block_size at the start, no need for extra length
	fwrite(out, 1, clen, out_fp);
	free(out);

	fastq_free(fq);
#endif
    }
    free(in);

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    if (hts_tpool_dispatch2(p, q, encode_thread, j, -1) != 0)
	goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
	enc_dec_job *j = hts_tpool_result_data(r);
	if (j->eof) {
	    end = 1;
	} else {
	    append_timings(t, &j->t, arg->verbose);
	    
	    // Track block offset for index
	    if (idx.nblocks >= idx_capacity) {
		idx_capacity *= 2;
		index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
		if (!new_entries) {
		    free(idx.entries);
		    hts_tpool_process_destroy(q);
		    hts_tpool_destroy(p);
		    return -1;
		}
		idx.entries = new_entries;
	    }
	    idx.entries[idx.nblocks].offset = ftell(out_fp);
	    idx.entries[idx.nblocks].usize = j->usize;
	    idx.entries[idx.nblocks].nrecords = j->nrecords;
	    idx.nblocks++;
	    
	    // Block now includes block_size at the start, no need for extra length
	    fwrite(j->comp, 1, j->clen, out_fp);
	    free(j->comp);
	}
	hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    // Write index at end of file
    uint64_t index_offset = ftell(out_fp);
    if (write_index(out_fp, &idx) < 0) {
        free(idx.entries);
        return -1;
    }
    
    // Update header with index offset
    fseek(out_fp, FQZ5_MAGIC_LEN, SEEK_SET);
    fwrite(&index_offset, 1, 8, out_fp);
    fseek(out_fp, 0, SEEK_END);
    
    free(idx.entries);
    return 0;

 err:
    free(idx.entries);
    return -1;
}

// Encode using gzFile (supports both plain and gzipped files)
int encode_gzip(gzFile in_fp, FILE *out_fp, fqz_gparams *gp, opts *arg,
           timings *t) {
    int rans_methods = (1<<RANS0) | (1<<RANS1) | (1<<RANS129) | (1<<RANS193);

    // Write file header with magic number
    if (write_header(out_fp) < 0)
        return -1;
    
    // Initialize index
    fqz5_index idx = {0};
    uint32_t idx_capacity = 1000;
    idx.entries = calloc(idx_capacity, sizeof(*idx.entries));
    if (!idx.entries)
        return -1;

    // Name
    if (arg->nauto) {
        method_avail[SEC_NAME] = arg->nauto;
    } else {
        if (arg->nstrat == 1)
            method_avail[SEC_NAME] |= 1<<(TOK3_3 + arg->nlevel/2-1);
        else if (arg->nstrat == 2)
            method_avail[SEC_NAME] |= 1<<(TOK3_3_LZP + arg->nlevel/2-1);
        else
            method_avail[SEC_NAME] = 1<<TLZP3;
    }
    
    // Seq
    if (arg->scustom) {
        method_avail[SEC_SEQ] = 1<<SEQ_CUSTOM;
    } else {
        if (arg->sauto)
            method_avail[SEC_SEQ] = arg->sauto;
        else if (arg->sstrat == 1)
            method_avail[SEC_SEQ] = 1<<SEQ_CUSTOM;

        if (!method_avail[SEC_SEQ])
            method_avail[SEC_SEQ] = rans_methods;
    }

    // Qual
    if (arg->qauto) {
        method_avail[SEC_QUAL] = arg->qauto;
    } else {
        if (arg->qstrat == 1) {
            if (arg->qlevel == 4)
                method_avail[SEC_QUAL] = FQZ4;
            else if (arg->qlevel == 3)
                method_avail[SEC_QUAL] = FQZ3;
            else if (arg->qlevel == 1)
                method_avail[SEC_QUAL] = FQZ1;
            else if (arg->qlevel == 2)
                method_avail[SEC_QUAL] = FQZ2;
            else
                method_avail[SEC_QUAL] = FQZ0;
        } else {
            method_avail[SEC_QUAL] = rans_methods;
        }
    }

#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n, 0);
    hts_tpool_result *r;
    enc_dec_job *j, *jr;
#endif

    int eof_flag = 0;

    while (!eof_flag) {
        fastq *fq = load_seqs_kseq(in_fp, arg->blk_size, &eof_flag);

        if (!fq) {
            free(idx.entries);
            return -1;
        }

        if (!fq->num_records) {
            fastq_free(fq);
            break;
        }

#ifdef THREADED
        // Dispatch a job
        j = calloc(1, sizeof(*j));
        j->gp = gp;
        j->arg = arg;
        memset(&j->t, 0, sizeof(j->t));
        j->fq = fq;
        j->eof = 0;
        j->usize = fq->seq_len;  // Total uncompressed bases
        j->nrecords = fq->num_records;

        int ret = -1;
        while (ret == -1) {
            // Always dispatch, going over-size on queue
            if ((ret=hts_tpool_dispatch2(p, q, encode_thread, j, -1)) != 0) {
                if (errno != EAGAIN)
                    goto err;
            }

            // Check for a result.
            // If input queue is oversize then do this blocking so we
            // don't grow input queue indefinitely.
            do {
                if (hts_tpool_dispatch_would_block(p, q)) {
                    r = hts_tpool_next_result_wait(q);
                } else {
                    r = hts_tpool_next_result(q);
                }
                if (r) {
                    jr = hts_tpool_result_data(r);
                    if (jr->eof) {
                        end = 1;
                    } else {
                        append_timings(t, &jr->t, arg->verbose);
                        
                        // Track block offset for index
                        if (idx.nblocks >= idx_capacity) {
                            idx_capacity *= 2;
                            index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
                            if (!new_entries) {
                                free(idx.entries);
                                goto err;
                            }
                            idx.entries = new_entries;
                        }
                        idx.entries[idx.nblocks].offset = ftell(out_fp);
                        idx.entries[idx.nblocks].usize = jr->usize;
                        idx.entries[idx.nblocks].nrecords = jr->nrecords;
                        idx.nblocks++;
                        
                        // Block now includes block_size at the start, no need for extra length
                        fwrite(jr->comp, 1, jr->clen, out_fp);
                        free(jr->comp);
                    }
                    hts_tpool_delete_result(r, 1);
                }
            } while (r && hts_tpool_dispatch_would_block(p, q));
        }
#else
        uint32_t clen;
        t->nblock++;
        
        // Track block offset for index
        if (idx.nblocks >= idx_capacity) {
            idx_capacity *= 2;
            index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
            if (!new_entries) {
                free(idx.entries);
                return -1;
            }
            idx.entries = new_entries;
        }
        idx.entries[idx.nblocks].offset = ftell(out_fp);
        idx.entries[idx.nblocks].usize = fq->seq_len;
        idx.entries[idx.nblocks].nrecords = fq->num_records;
        idx.nblocks++;
        
        char *out = encode_block(gp, arg, fq, t, &clen);
        // Block now includes block_size at the start, no need for extra length
        fwrite(out, 1, clen, out_fp);
        free(out);

        fastq_free(fq);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    if (hts_tpool_dispatch2(p, q, encode_thread, j, -1) != 0)
        goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
        enc_dec_job *jr = hts_tpool_result_data(r);
        if (jr->eof) {
            end = 1;
        } else {
            append_timings(t, &jr->t, arg->verbose);
            
            // Track block offset for index
            if (idx.nblocks >= idx_capacity) {
                idx_capacity *= 2;
                index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
                if (!new_entries) {
                    free(idx.entries);
                    hts_tpool_process_destroy(q);
                    hts_tpool_destroy(p);
                    return -1;
                }
                idx.entries = new_entries;
            }
            idx.entries[idx.nblocks].offset = ftell(out_fp);
            idx.entries[idx.nblocks].usize = jr->usize;
            idx.entries[idx.nblocks].nrecords = jr->nrecords;
            idx.nblocks++;
            
            // Block now includes block_size at the start, no need for extra length
            fwrite(jr->comp, 1, jr->clen, out_fp);
            free(jr->comp);
        }
        hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    // Write index at end of file
    uint64_t index_offset = ftell(out_fp);
    if (write_index(out_fp, &idx) < 0) {
        free(idx.entries);
        return -1;
    }
    
    // Update header with index offset
    fseek(out_fp, FQZ5_MAGIC_LEN, SEEK_SET);
    fwrite(&index_offset, 1, 8, out_fp);
    fseek(out_fp, 0, SEEK_END);
    
    free(idx.entries);
    return 0;

 err:
    free(idx.entries);
    return -1;
}

// Encode two paired FASTQ files with interleaving
int encode_interleaved(gzFile in_fp1, gzFile in_fp2, FILE *out_fp, fqz_gparams *gp, opts *arg,
           timings *t) {
    int rans_methods = (1<<RANS0) | (1<<RANS1) | (1<<RANS129) | (1<<RANS193);

    // Write file header with magic number
    if (write_header(out_fp) < 0)
        return -1;
    
    // Initialize index
    fqz5_index idx = {0};
    uint32_t idx_capacity = 1000;
    idx.entries = calloc(idx_capacity, sizeof(*idx.entries));
    if (!idx.entries)
        return -1;

    // Name
    if (arg->nauto) {
        method_avail[SEC_NAME] = arg->nauto;
    } else {
        if (arg->nstrat == 1)
            method_avail[SEC_NAME] |= 1<<(TOK3_3 + arg->nlevel/2-1);
        else if (arg->nstrat == 2)
            method_avail[SEC_NAME] |= 1<<(TOK3_3_LZP + arg->nlevel/2-1);
        else
            method_avail[SEC_NAME] = 1<<TLZP3;
    }
    
    // Seq
    if (arg->scustom) {
        method_avail[SEC_SEQ] = 1<<SEQ_CUSTOM;
    } else {
        if (arg->sauto)
            method_avail[SEC_SEQ] = arg->sauto;
        else if (arg->sstrat == 1)
            method_avail[SEC_SEQ] = 1<<SEQ_CUSTOM;

        if (!method_avail[SEC_SEQ])
            method_avail[SEC_SEQ] = rans_methods;
    }

    // Qual
    if (arg->qauto) {
        method_avail[SEC_QUAL] = arg->qauto;
    } else {
        if (arg->qstrat == 1) {
            if (arg->qlevel == 4)
                method_avail[SEC_QUAL] = FQZ4;
            else if (arg->qlevel == 3)
                method_avail[SEC_QUAL] = FQZ3;
            else if (arg->qlevel == 1)
                method_avail[SEC_QUAL] = FQZ1;
            else if (arg->qlevel == 2)
                method_avail[SEC_QUAL] = FQZ2;
            else
                method_avail[SEC_QUAL] = FQZ0;
        } else {
            method_avail[SEC_QUAL] = rans_methods;
        }
    }

#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n, 0);
    hts_tpool_result *r;
    enc_dec_job *j, *jr;
#endif

    int eof_flag = 0;

    while (!eof_flag) {
        fastq *fq = load_seqs_interleaved(in_fp1, in_fp2, arg->blk_size, &eof_flag);

        if (!fq) {
            free(idx.entries);
            return -1;
        }

        if (!fq->num_records) {
            fastq_free(fq);
            break;
        }

#ifdef THREADED
        // Dispatch a job
        j = calloc(1, sizeof(*j));
        j->gp = gp;
        j->arg = arg;
        memset(&j->t, 0, sizeof(j->t));
        j->fq = fq;
        j->eof = 0;
        j->usize = fq->seq_len;  // Total uncompressed bases
        j->nrecords = fq->num_records;

        int ret = -1;
        while (ret == -1) {
            // Always dispatch, going over-size on queue
            if ((ret=hts_tpool_dispatch2(p, q, encode_thread, j, -1)) != 0) {
                if (errno != EAGAIN)
                    goto err;
            }

            // Check for a result.
            // If input queue is oversize then do this blocking so we
            // don't grow input queue indefinitely.
            do {
                if (hts_tpool_dispatch_would_block(p, q)) {
                    r = hts_tpool_next_result_wait(q);
                } else {
                    r = hts_tpool_next_result(q);
                }
                if (r) {
                    jr = hts_tpool_result_data(r);
                    if (jr->eof) {
                        end = 1;
                    } else {
                        append_timings(t, &jr->t, arg->verbose);
                        
                        // Track block offset for index
                        if (idx.nblocks >= idx_capacity) {
                            idx_capacity *= 2;
                            index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
                            if (!new_entries) {
                                free(idx.entries);
                                goto err;
                            }
                            idx.entries = new_entries;
                        }
                        idx.entries[idx.nblocks].offset = ftell(out_fp);
                        idx.entries[idx.nblocks].usize = jr->usize;
                        idx.entries[idx.nblocks].nrecords = jr->nrecords;
                        idx.nblocks++;
                        
                        // Block now includes block_size at the start, no need for extra length
                        fwrite(jr->comp, 1, jr->clen, out_fp);
                        free(jr->comp);
                    }
                    hts_tpool_delete_result(r, 1);
                }
            } while (r && hts_tpool_dispatch_would_block(p, q));
        }
#else
        uint32_t clen;
        t->nblock++;
        
        // Track block offset for index
        if (idx.nblocks >= idx_capacity) {
            idx_capacity *= 2;
            index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
            if (!new_entries) {
                free(idx.entries);
                return -1;
            }
            idx.entries = new_entries;
        }
        idx.entries[idx.nblocks].offset = ftell(out_fp);
        idx.entries[idx.nblocks].usize = fq->seq_len;
        idx.entries[idx.nblocks].nrecords = fq->num_records;
        idx.nblocks++;
        
        char *out = encode_block(gp, arg, fq, t, &clen);
        // Block now includes block_size at the start, no need for extra length
        fwrite(out, 1, clen, out_fp);
        free(out);

        fastq_free(fq);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    if (hts_tpool_dispatch2(p, q, encode_thread, j, -1) != 0)
        goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
        enc_dec_job *jr = hts_tpool_result_data(r);
        if (jr->eof) {
            end = 1;
        } else {
            append_timings(t, &jr->t, arg->verbose);
            
            // Track block offset for index
            if (idx.nblocks >= idx_capacity) {
                idx_capacity *= 2;
                index_entry *new_entries = realloc(idx.entries, idx_capacity * sizeof(*idx.entries));
                if (!new_entries) {
                    free(idx.entries);
                    hts_tpool_process_destroy(q);
                    hts_tpool_destroy(p);
                    return -1;
                }
                idx.entries = new_entries;
            }
            idx.entries[idx.nblocks].offset = ftell(out_fp);
            idx.entries[idx.nblocks].usize = jr->usize;
            idx.entries[idx.nblocks].nrecords = jr->nrecords;
            idx.nblocks++;
            
            // Block now includes block_size at the start, no need for extra length
            fwrite(jr->comp, 1, jr->clen, out_fp);
            free(jr->comp);
        }
        hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    // Write index at end of file
    uint64_t index_offset = ftell(out_fp);
    if (write_index(out_fp, &idx) < 0) {
        free(idx.entries);
        return -1;
    }
    
    // Update header with index offset
    fseek(out_fp, FQZ5_MAGIC_LEN, SEEK_SET);
    fwrite(&index_offset, 1, 8, out_fp);
    fseek(out_fp, 0, SEEK_END);
    
    free(idx.entries);
    return 0;

 err:
    free(idx.entries);
    return -1;
}

int output_fastq(FILE *out_fp, fastq *fq, int plus_name) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;
    char *qp = fq->qual_buf;

#if 1
    // A bit faster sometimes.
    // Calculate buffer size - if plus_name, need space for name on third line
    // Base: name_len (with nulls) + seq_len + qual_len + num_records*5 (for @,\n,\n,+,\n,\n per record)
    // With plus_name: add name_len (names duplicated on line 3) + num_records (for safety)
    int len = fq->name_len + fq->seq_len + fq->qual_len + fq->num_records*5;
    if (plus_name)
        len += fq->name_len + fq->num_records;
    char *buf = malloc(len), *cp = buf;

    for (int i = 0; i < fq->num_records; i++) {
	char *name_start = np;  // save start of name for third line
	*cp++ = '@';
	while ((*cp++ = *np++))
	    ;
	*--cp = '\n'; cp++;
	memcpy(cp, sp, fq->len[i]);
	cp += fq->len[i];
	sp += fq->len[i];
	*cp++ = '\n';
	*cp++ = '+';
	if (plus_name) {
	    // Copy name to third line
	    char *tmp = name_start;
	    while ((*cp++ = *tmp++))
		;
	    *--cp = '\n'; cp++;
	} else {
	    *cp++ = '\n';
	}
	memcpy(cp, qp, fq->len[i]);
	cp += fq->len[i];
	qp += fq->len[i];
	*cp++ = '\n';
    }
    fwrite(buf, 1, cp-buf, out_fp);
    free(buf);
#else
    for (int i = 0; i < fq->num_records; i++) {
	if (plus_name) {
	    char *name_copy = np;
	    fprintf(out_fp, "@%s\n%.*s\n+%s\n%.*s\n",
		    np, fq->len[i], sp, name_copy, fq->len[i], qp);
	} else {
	    fprintf(out_fp, "@%s\n%.*s\n+\n%.*s\n",
		    np, fq->len[i], sp, fq->len[i], qp);
	}
	np += strlen(np)+1;
	sp += fq->len[i];
	qp += fq->len[i];
    }
#endif

    return 0;
}

// Output FASTA format
int output_fasta(FILE *out_fp, fastq *fq) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;


    // Use fprintf for safety
    for (int i = 0; i < fq->num_records; i++) {
	fprintf(out_fp, ">%s\n", np);
	fwrite(sp, 1, fq->len[i], out_fp);
	fputc('\n', out_fp);
	np += strlen(np) + 1;
	sp += fq->len[i];
    }
    return 0;
}

int output_fasta_gzip(gzFile out_fp, fastq *fq) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;

    // Use gzprintf for safety
    for (int i = 0; i < fq->num_records; i++) {
	gzprintf(out_fp, ">%s\n", np);
	gzwrite(out_fp, sp, fq->len[i]);
	gzputc(out_fp, '\n');
	np += strlen(np) + 1;
	sp += fq->len[i];
    }
    return 0;
}

// Output FASTA format to two files (deinterleaved)
int output_fasta_deinterleaved(FILE *out_fp1, FILE *out_fp2, fastq *fq) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;

    for (int i = 0; i < fq->num_records; i++) {
	FILE *out_fp = (i & 1) ? out_fp2 : out_fp1;
	fprintf(out_fp, ">%s\n", np);
	fwrite(sp, 1, fq->len[i], out_fp);
	fputc('\n', out_fp);
	np += strlen(np) + 1;
	sp += fq->len[i];
    }
    return 0;
}

// Output FASTA format to two gzip files (deinterleaved)
int output_fasta_gzip_deinterleaved(gzFile out_fp1, gzFile out_fp2, fastq *fq) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;

    for (int i = 0; i < fq->num_records; i++) {
	gzFile out_fp = (i & 1) ? out_fp2 : out_fp1;
	gzprintf(out_fp, ">%s\n", np);
	gzwrite(out_fp, sp, fq->len[i]);
	gzputc(out_fp, '\n');
	np += strlen(np) + 1;
	sp += fq->len[i];
    }
    return 0;
}

int output_fastq_gzip(gzFile out_fp, fastq *fq, int plus_name) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;
    char *qp = fq->qual_buf;

    // Build buffer and write at once
    // Calculate buffer size - if plus_name, need space for name on third line
    // Base: name_len (with nulls) + seq_len + qual_len + num_records*5 (for @,\n,\n,+,\n,\n per record)
    // With plus_name: add name_len (names duplicated on line 3) + num_records (for safety)
    int len = fq->name_len + fq->seq_len + fq->qual_len + fq->num_records*5;
    if (plus_name)
        len += fq->name_len + fq->num_records;
    char *buf = malloc(len), *cp = buf;

    for (int i = 0; i < fq->num_records; i++) {
        char *name_start = np;  // save start of name for third line
        *cp++ = '@';
        while ((*cp++ = *np++))
            ;
        *--cp = '\n'; cp++;
        memcpy(cp, sp, fq->len[i]);
        cp += fq->len[i];
        sp += fq->len[i];
        *cp++ = '\n';
        *cp++ = '+';
        if (plus_name) {
            // Copy name to third line
            char *tmp = name_start;
            while ((*cp++ = *tmp++))
                ;
            *--cp = '\n'; cp++;
        } else {
            *cp++ = '\n';
        }
        memcpy(cp, qp, fq->len[i]);
        cp += fq->len[i];
        qp += fq->len[i];
        *cp++ = '\n';
    }
    gzwrite(out_fp, buf, cp-buf);
    free(buf);

    return 0;
}

// Output interleaved FASTQ to two separate files (deinterleaving)
int output_fastq_deinterleaved(FILE *out_fp1, FILE *out_fp2, fastq *fq, int plus_name) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;
    char *qp = fq->qual_buf;

    // Build separate buffers for R1 and R2
    // Calculate buffer size - if plus_name, need space for name on third line
    // Base: name_len (with nulls) + seq_len + qual_len + num_records*5 (for @,\n,\n,+,\n,\n per record)
    // With plus_name: add name_len (names duplicated on line 3) + num_records (for safety)
    int len = fq->name_len + fq->seq_len + fq->qual_len + fq->num_records*5;
    if (plus_name)
        len += fq->name_len + fq->num_records;
    char *buf1 = malloc(len);
    char *buf2 = malloc(len);
    if (!buf1 || !buf2) {
        free(buf1);
        free(buf2);
        fprintf(stderr, "Failed to allocate memory for deinterleaving\n");
        return -1;
    }
    char *cp1 = buf1;
    char *cp2 = buf2;

    for (int i = 0; i < fq->num_records; i++) {
        // Determine which buffer to use based on position
        // Even indices (0, 2, 4, ...) go to file 1 (R1)
        // Odd indices (1, 3, 5, ...) go to file 2 (R2)
        char **cpp = (i % 2 == 1) ? &cp2 : &cp1;
        char *cp = *cpp;
        
        char *name_start = np;  // save start of name for third line
        *cp++ = '@';
        while ((*cp++ = *np++))
            ;
        *--cp = '\n'; cp++;
        memcpy(cp, sp, fq->len[i]);
        cp += fq->len[i];
        sp += fq->len[i];
        *cp++ = '\n';
        *cp++ = '+';
        if (plus_name) {
            // Copy name to third line
            char *tmp = name_start;
            while ((*cp++ = *tmp++))
                ;
            *--cp = '\n'; cp++;
        } else {
            *cp++ = '\n';
        }
        memcpy(cp, qp, fq->len[i]);
        cp += fq->len[i];
        qp += fq->len[i];
        *cp++ = '\n';
        
        *cpp = cp;
    }
    
    fwrite(buf1, 1, cp1-buf1, out_fp1);
    fwrite(buf2, 1, cp2-buf2, out_fp2);
    free(buf1);
    free(buf2);

    return 0;
}

// Output interleaved FASTQ to two separate gzipped files (deinterleaving)
int output_fastq_gzip_deinterleaved(gzFile out_fp1, gzFile out_fp2, fastq *fq, int plus_name) {
    char *np = fq->name_buf;
    char *sp = fq->seq_buf;
    char *qp = fq->qual_buf;

    // Build separate buffers for R1 and R2
    // Calculate buffer size - if plus_name, need space for name on third line
    // Base: name_len (with nulls) + seq_len + qual_len + num_records*5 (for @,\n,\n,+,\n,\n per record)
    // With plus_name: add name_len (names duplicated on line 3) + num_records (for safety)
    int len = fq->name_len + fq->seq_len + fq->qual_len + fq->num_records*5;
    if (plus_name)
        len += fq->name_len + fq->num_records;
    char *buf1 = malloc(len);
    char *buf2 = malloc(len);
    if (!buf1 || !buf2) {
        free(buf1);
        free(buf2);
        fprintf(stderr, "Failed to allocate memory for deinterleaving\n");
        return -1;
    }
    char *cp1 = buf1;
    char *cp2 = buf2;

    for (int i = 0; i < fq->num_records; i++) {
        // Determine which buffer to use based on position
        // Even indices (0, 2, 4, ...) go to file 1 (R1)
        // Odd indices (1, 3, 5, ...) go to file 2 (R2)
        char **cpp = (i % 2 == 1) ? &cp2 : &cp1;
        char *cp = *cpp;
        
        char *name_start = np;  // save start of name for third line
        *cp++ = '@';
        while ((*cp++ = *np++))
            ;
        *--cp = '\n'; cp++;
        memcpy(cp, sp, fq->len[i]);
        cp += fq->len[i];
        sp += fq->len[i];
        *cp++ = '\n';
        *cp++ = '+';
        if (plus_name) {
            // Copy name to third line
            char *tmp = name_start;
            while ((*cp++ = *tmp++))
                ;
            *--cp = '\n'; cp++;
        } else {
            *cp++ = '\n';
        }
        memcpy(cp, qp, fq->len[i]);
        cp += fq->len[i];
        qp += fq->len[i];
        *cp++ = '\n';
        
        *cpp = cp;
    }
    
    gzwrite(out_fp1, buf1, cp1-buf1);
    gzwrite(out_fp2, buf2, cp2-buf2);
    free(buf1);
    free(buf2);

    return 0;
}

static void *decode_thread(void *arg) {
    enc_dec_job *j = (enc_dec_job *)arg;
    if (j->eof)
	return j;

    j->fq = decode_block(j->comp, j->clen, &j->t, j->file_version);
    free(j->comp);
    return j;
}

int decode(FILE *in_fp, FILE *out_fp, opts *arg, timings *t) {
    uint64_t index_offset;
    int header_result = read_header(in_fp, &index_offset);
    if (header_result < 0)
        return -1;
    // header_result: 0 = v1.1 (with CRC), 1 = v1.0 (no CRC), 2 = old format
    int file_version = header_result;
    
#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n, 0);
    hts_tpool_result *r;
    enc_dec_job *j;
#endif

    for (;;) {
	// Check if we've reached the index
	uint64_t current_pos = ftell(in_fp);
	if (index_offset > 0 && current_pos >= index_offset)
	    break;
	    
	// Load next compressed block
	int i;
	uint32_t block_size;
	unsigned char *comp;

	// Read block size (first 4 bytes of block)
	if (fread(&block_size, 1, 4, in_fp) != 4)
	    break;

	// Allocate buffer for entire block (including the block_size field we just read)
	uint32_t c_len = block_size + 4;  // +4 for the block_size field itself
	comp = malloc(c_len);
	if (!comp)
	    return -1;
	
	// Store the block_size we already read at the start
	*(uint32_t *)comp = block_size;
	
	// Read the rest of the block (block_size bytes)
	if (fread(comp + 4, 1, block_size, in_fp) != block_size) {
	    free(comp);
	    return -1;
	}

#ifdef THREADED
	// Dispatch a job
	j = calloc(1, sizeof(*j));
	memset(&j->t, 0, sizeof(j->t));
	j->comp = comp;
	j->clen = c_len;
	j->fq = NULL;
	j->eof = 0;
	j->file_version = file_version;

	// Always put on queue, even if over queue size
	if (hts_tpool_dispatch2(p, q, decode_thread, j, -1) != 0)
	    goto err;

	// Check for a result.
	do {
	    if (hts_tpool_dispatch_would_block(p, q)) {
		//fprintf(stderr, "would block: in=%d\n", q->n_input);
		r = hts_tpool_next_result_wait(q);
	    } else {
		//fprintf(stderr, "would pass : in=%d\n", q->n_input);
		r = hts_tpool_next_result(q);
	    }
	    if (r) {
		j = hts_tpool_result_data(r);
		if (j->eof) {
		    end = 1;
		} else {
		    if (!j->fq) {
			fprintf(stderr, "ERROR: Failed to decode block\n");
			hts_tpool_delete_result(r, 1);
			goto err;
		    }
		    append_timings(t, &j->t, arg->verbose);
		    if (j->fq->is_fasta)
			output_fasta(out_fp, j->fq);
		    else
			output_fastq(out_fp, j->fq, arg->plus_name);
		    fastq_free(j->fq);
		}
		hts_tpool_delete_result(r, 1);
	    }
	} while (r && hts_tpool_dispatch_would_block(p, q));
#else
	t->nblock++;
	fastq *fq = decode_block(comp, c_len, t, file_version);
	if (!fq) {
	    fprintf(stderr, "ERROR: Failed to decode block\n");
	    free(comp);
	    return -1;
	}

	// ----------
	// Convert back to fastq
	char *np = fq->name_buf;
	char *sp = fq->seq_buf;
	char *qp = fq->qual_buf;

	for (i = 0; i < fq->num_records; i++) {
	    fprintf(out_fp, "@%s\n%.*s\n+\n%.*s\n",
		    np, fq->len[i], sp, fq->len[i], qp);
	    np += strlen(np)+1;
	    sp += fq->len[i];
	    qp += fq->len[i];
	}

	fastq_free(fq);
	free(comp);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    j->file_version = file_version;
    j->file_version = file_version;
    j->eof = 1;
    if (hts_tpool_dispatch(p, q, decode_thread, j) != 0)
	goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
	enc_dec_job *j = hts_tpool_result_data(r);
	if (j->eof) {
	    end = 1;
	} else {
	    if (!j->fq) {
		fprintf(stderr, "ERROR: Failed to decode block\n");
		hts_tpool_delete_result(r, 1);
		goto err;
	    }
	    append_timings(t, &j->t, arg->verbose);
	    if (j->fq->is_fasta)
		output_fasta(out_fp, j->fq);
	    else
		output_fastq(out_fp, j->fq, arg->plus_name);
	    fastq_free(j->fq);
	}
	hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    return 0;

 err:
    // FIXME: tidy up
    return -1;
}

int decode_gzip(FILE *in_fp, gzFile out_fp, opts *arg, timings *t) {
    uint64_t index_offset;
    int header_result = read_header(in_fp, &index_offset);
    if (header_result < 0)
        return -1;
    // header_result: 0 = v1.1 (with CRC), 1 = v1.0 (no CRC), 2 = old format
    int file_version = header_result;
    
#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n, 0);
    hts_tpool_result *r;
    enc_dec_job *j;
#endif

    for (;;) {
        // Check if we've reached the index
        uint64_t current_pos = ftell(in_fp);
        if (index_offset > 0 && current_pos >= index_offset)
            break;
            
        // Load next compressed block
        int i;
        uint32_t block_size;
        unsigned char *comp;

        // Read block size (first 4 bytes of block)
        if (fread(&block_size, 1, 4, in_fp) != 4)
            break;

        // Allocate buffer for entire block (including the block_size field we just read)
        uint32_t c_len = block_size + 4;  // +4 for the block_size field itself
        comp = malloc(c_len);
        if (!comp)
            return -1;
        
        // Store the block_size we already read at the start
        *(uint32_t *)comp = block_size;
        
        // Read the rest of the block (block_size bytes)
        if (fread(comp + 4, 1, block_size, in_fp) != block_size) {
            free(comp);
            return -1;
        }

#ifdef THREADED
        // Dispatch a job
        j = calloc(1, sizeof(*j));
        memset(&j->t, 0, sizeof(j->t));
        j->comp = comp;
        j->clen = c_len;
        j->fq = NULL;
        j->eof = 0;
        j->file_version = file_version;

        // Always put on queue, even if over queue size
        if (hts_tpool_dispatch2(p, q, decode_thread, j, -1) != 0)
            goto err;

        // Check for a result.
        do {
            if (hts_tpool_dispatch_would_block(p, q)) {
                r = hts_tpool_next_result_wait(q);
            } else {
                r = hts_tpool_next_result(q);
            }
            if (r) {
                j = hts_tpool_result_data(r);
                if (j->eof) {
                    end = 1;
                } else {
                    if (!j->fq) {
                        fprintf(stderr, "ERROR: Failed to decode block\n");
                        hts_tpool_delete_result(r, 1);
                        goto err;
                    }
                    append_timings(t, &j->t, arg->verbose);
                    if (j->fq->is_fasta)
                        output_fasta_gzip(out_fp, j->fq);
                    else
                        output_fastq_gzip(out_fp, j->fq, arg->plus_name);
                    fastq_free(j->fq);
                }
                hts_tpool_delete_result(r, 1);
            }
        } while (r && hts_tpool_dispatch_would_block(p, q));
#else
        t->nblock++;
        fastq *fq = decode_block(comp, c_len, t, file_version);
        if (!fq) {
            fprintf(stderr, "ERROR: Failed to decode block\n");
            free(comp);
            return -1;
        }

        output_fastq_gzip(out_fp, fq, arg->plus_name);

        fastq_free(fq);
        free(comp);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    j->file_version = file_version;
    if (hts_tpool_dispatch(p, q, decode_thread, j) != 0)
        goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
        enc_dec_job *j = hts_tpool_result_data(r);
        if (j->eof) {
            end = 1;
        } else {
            if (!j->fq) {
                fprintf(stderr, "ERROR: Failed to decode block\n");
                hts_tpool_delete_result(r, 1);
                goto err;
            }
            append_timings(t, &j->t, arg->verbose);
            output_fastq_gzip(out_fp, j->fq, arg->plus_name);
            fastq_free(j->fq);
        }
        hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    return 0;

 err:
    // FIXME: tidy up
    return -1;
}

// Decode with deinterleaving to two separate files
int decode_deinterleaved(FILE *in_fp, FILE *out_fp1, FILE *out_fp2, opts *arg, timings *t) {
    uint64_t index_offset;
    int header_result = read_header(in_fp, &index_offset);
    if (header_result < 0)
        return -1;
    // header_result: 0 = v1.1 (with CRC), 1 = v1.0 (no CRC), 2 = old format
    int file_version = header_result;
    
#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n, 0);
    hts_tpool_result *r;
    enc_dec_job *j;
#endif

    for (;;) {
        // Check if we've reached the index
        uint64_t current_pos = ftell(in_fp);
        if (index_offset > 0 && current_pos >= index_offset)
            break;
            
        // Load next compressed block
        int i;
        uint32_t block_size;
        unsigned char *comp;

        // Read block size (first 4 bytes of block)
        if (fread(&block_size, 1, 4, in_fp) != 4)
            break;

        // Allocate buffer for entire block (including the block_size field we just read)
        uint32_t c_len = block_size + 4;  // +4 for the block_size field itself
        comp = malloc(c_len);
        if (!comp)
            return -1;
        
        // Store the block_size we already read at the start
        *(uint32_t *)comp = block_size;
        
        // Read the rest of the block (block_size bytes)
        if (fread(comp + 4, 1, block_size, in_fp) != block_size) {
            free(comp);
            return -1;
        }

#ifdef THREADED
        // Dispatch a job
        j = calloc(1, sizeof(*j));
        memset(&j->t, 0, sizeof(j->t));
        j->comp = comp;
        j->clen = c_len;
        j->fq = NULL;
        j->eof = 0;
        j->file_version = file_version;

        // Always put on queue, even if over queue size
        if (hts_tpool_dispatch2(p, q, decode_thread, j, -1) != 0)
            goto err;

        // Check for a result.
        do {
            if (hts_tpool_dispatch_would_block(p, q)) {
                r = hts_tpool_next_result_wait(q);
            } else {
                r = hts_tpool_next_result(q);
            }
            if (r) {
                j = hts_tpool_result_data(r);
                if (j->eof) {
                    end = 1;
                } else {
                    if (!j->fq) {
                        fprintf(stderr, "ERROR: Failed to decode block\n");
                        hts_tpool_delete_result(r, 1);
                        goto err;
                    }
                    append_timings(t, &j->t, arg->verbose);
                    if (j->fq->is_fasta)
                        output_fasta_deinterleaved(out_fp1, out_fp2, j->fq);
                    else
                        if (j->fq->is_fasta)
                output_fasta_deinterleaved(out_fp1, out_fp2, j->fq);
            else
                output_fastq_deinterleaved(out_fp1, out_fp2, j->fq, arg->plus_name);
                    fastq_free(j->fq);
                }
                hts_tpool_delete_result(r, 1);
            }
        } while (r && hts_tpool_dispatch_would_block(p, q));
#else
        t->nblock++;
        fastq *fq = decode_block(comp, c_len, t, file_version);
        if (!fq) {
            fprintf(stderr, "ERROR: Failed to decode block\n");
            free(comp);
            return -1;
        }

        if (fq->is_fasta)
            output_fasta_deinterleaved(out_fp1, out_fp2, fq);
        else
            output_fastq_deinterleaved(out_fp1, out_fp2, fq, arg->plus_name);

        free(comp);
        fastq_free(fq);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    j->file_version = file_version;
    if (hts_tpool_dispatch2(p, q, decode_thread, j, -1) != 0)
        goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
        enc_dec_job *j = hts_tpool_result_data(r);
        if (j->eof) {
            end = 1;
        } else {
            if (!j->fq) {
                fprintf(stderr, "ERROR: Failed to decode block\n");
                hts_tpool_delete_result(r, 1);
                goto err;
            }
            append_timings(t, &j->t, arg->verbose);
            if (j->fq->is_fasta)
                output_fasta_deinterleaved(out_fp1, out_fp2, j->fq);
            else
                output_fastq_deinterleaved(out_fp1, out_fp2, j->fq, arg->plus_name);
            fastq_free(j->fq);
        }
        hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    return 0;

 err:
    // FIXME: tidy up
    return -1;
}

// Decode with deinterleaving to two separate gzipped files
int decode_gzip_deinterleaved(FILE *in_fp, gzFile out_fp1, gzFile out_fp2, opts *arg, timings *t) {
    uint64_t index_offset;
    int header_result = read_header(in_fp, &index_offset);
    if (header_result < 0)
        return -1;
    // header_result: 0 = v1.1 (with CRC), 1 = v1.0 (no CRC), 2 = old format
    int file_version = header_result;
    
#ifdef THREADED
    int n = arg->nthread, end = 0;
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n, 0);
    hts_tpool_result *r;
    enc_dec_job *j;
#endif

    for (;;) {
        // Check if we've reached the index
        uint64_t current_pos = ftell(in_fp);
        if (index_offset > 0 && current_pos >= index_offset)
            break;
            
        // Load next compressed block
        int i;
        uint32_t block_size;
        unsigned char *comp;

        // Read block size (first 4 bytes of block)
        if (fread(&block_size, 1, 4, in_fp) != 4)
            break;

        // Allocate buffer for entire block (including the block_size field we just read)
        uint32_t c_len = block_size + 4;  // +4 for the block_size field itself
        comp = malloc(c_len);
        if (!comp)
            return -1;
        
        // Store the block_size we already read at the start
        *(uint32_t *)comp = block_size;
        
        // Read the rest of the block (block_size bytes)
        if (fread(comp + 4, 1, block_size, in_fp) != block_size) {
            free(comp);
            return -1;
        }

#ifdef THREADED
        // Dispatch a job
        j = calloc(1, sizeof(*j));
        memset(&j->t, 0, sizeof(j->t));
        j->comp = comp;
        j->clen = c_len;
        j->fq = NULL;
        j->eof = 0;
        j->file_version = file_version;

        // Always put on queue, even if over queue size
        if (hts_tpool_dispatch2(p, q, decode_thread, j, -1) != 0)
            goto err;

        // Check for a result.
        do {
            if (hts_tpool_dispatch_would_block(p, q)) {
                r = hts_tpool_next_result_wait(q);
            } else {
                r = hts_tpool_next_result(q);
            }
            if (r) {
                j = hts_tpool_result_data(r);
                if (j->eof) {
                    end = 1;
                } else {
                    if (!j->fq) {
                        fprintf(stderr, "ERROR: Failed to decode block\n");
                        hts_tpool_delete_result(r, 1);
                        goto err;
                    }
                    append_timings(t, &j->t, arg->verbose);
                    if (j->fq->is_fasta)
                        output_fasta_gzip_deinterleaved(out_fp1, out_fp2, j->fq);
                    else
                        if (j->fq->is_fasta)
                output_fasta_gzip_deinterleaved(out_fp1, out_fp2, j->fq);
            else
                output_fastq_gzip_deinterleaved(out_fp1, out_fp2, j->fq, arg->plus_name);
                    fastq_free(j->fq);
                }
                hts_tpool_delete_result(r, 1);
            }
        } while (r && hts_tpool_dispatch_would_block(p, q));
#else
        t->nblock++;
        fastq *fq = decode_block(comp, c_len, t, file_version);
        if (!fq) {
            fprintf(stderr, "ERROR: Failed to decode block\n");
            free(comp);
            return -1;
        }

        if (fq->is_fasta)
            output_fasta_gzip_deinterleaved(out_fp1, out_fp2, fq);
        else
            output_fastq_gzip_deinterleaved(out_fp1, out_fp2, fq, arg->plus_name);

        free(comp);
        fastq_free(fq);
#endif
    }

#ifdef THREADED
    j = malloc(sizeof(*j));
    j->eof = 1;
    j->file_version = file_version;
    if (hts_tpool_dispatch2(p, q, decode_thread, j, -1) != 0)
        goto err;

    // End of input, so work through remaining results
    while (!end && (r = hts_tpool_next_result_wait(q))) {
        enc_dec_job *j = hts_tpool_result_data(r);
        if (j->eof) {
            end = 1;
        } else {
            if (!j->fq) {
                fprintf(stderr, "ERROR: Failed to decode block\n");
                hts_tpool_delete_result(r, 1);
                goto err;
            }
            append_timings(t, &j->t, arg->verbose);
            if (j->fq->is_fasta)
                output_fasta_gzip_deinterleaved(out_fp1, out_fp2, j->fq);
            else
                output_fastq_gzip_deinterleaved(out_fp1, out_fp2, j->fq, arg->plus_name);
            fastq_free(j->fq);
        }
        hts_tpool_delete_result(r, 1);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
#endif

    return 0;

 err:
    // FIXME: tidy up
    return -1;
}

// Inspect FQZ5 file and display comprehensive information
int inspect_file(FILE *in_fp, opts *arg) {
    uint64_t index_offset;
    int header_result = read_header(in_fp, &index_offset);
    
    // Determine file version
    const char *version_str;
    int has_crc = 0;
    if (header_result < 0) {
        fprintf(stderr, "Error: Failed to read file header\n");
        return -1;
    } else if (header_result == 0) {
        version_str = "1.1 (current)";
        has_crc = 1;
    } else if (header_result == 1) {
        version_str = "1.0 (legacy)";
        has_crc = 0;
    } else {
        version_str = "pre-1.0 (legacy, no header)";
        has_crc = 0;
    }
    
    printf("FQZ5 File Inspection\n");
    printf("====================\n\n");
    
    // Display version
    printf("Format Version:      %s\n", version_str);
    
    // Get file size
    fseek(in_fp, 0, SEEK_END);
    uint64_t file_size = ftell(in_fp);
    fseek(in_fp, FQZ5_MAGIC_LEN + 8, SEEK_SET); // Skip header
    if (header_result == 2) {
        fseek(in_fp, 0, SEEK_SET); // Old format, no header
    }
    
    printf("Compressed Size:     %llu bytes (%.2f MB)\n", 
           (unsigned long long)file_size, file_size / 1048576.0);
    
    // Parse blocks to gather statistics
    int nblocks = 0;
    uint64_t total_uncompressed = 0;
    uint64_t total_records = 0;
    uint64_t total_compressed_data = 0;
    int integrity_errors = 0;
    
    for (;;) {
        // Check if we've reached the index
        uint64_t current_pos = ftell(in_fp);
        if (index_offset > 0 && current_pos >= index_offset)
            break;
            
        // Read block size (first 4 bytes of block)
        uint32_t block_size;
        if (fread(&block_size, 1, 4, in_fp) != 4)
            break;
        
        // Validate block size to prevent underflow
        if (has_crc && block_size < 8) {
            fprintf(stderr, "Warning: Invalid block size %u in block %d (too small)\n", block_size, nblocks);
            break;
        } else if (!has_crc && block_size < 4) {
            fprintf(stderr, "Warning: Invalid block size %u in block %d (too small)\n", block_size, nblocks);
            break;
        }
        
        // Read num_records
        uint32_t num_records;
        if (fread(&num_records, 1, 4, in_fp) != 4) {
            fprintf(stderr, "Warning: Failed to read num_records in block %d\n", nblocks);
            break;
        }
        
        total_records += num_records;
        total_compressed_data += block_size;
        
        // If we have CRC, read and verify it
        if (has_crc) {
            uint32_t stored_crc;
            if (fread(&stored_crc, 1, 4, in_fp) != 4) {
                fprintf(stderr, "Warning: Failed to read CRC in block %d\n", nblocks);
                break;
            }
            
            // Read the rest of the block data
            uint32_t data_size = block_size - 8;  // Subtract num_records (4) and CRC (4)
            unsigned char *data = malloc(data_size);
            if (!data) {
                fprintf(stderr, "Warning: Failed to allocate memory for block %d\n", nblocks);
                break;
            }
            
            if (fread(data, 1, data_size, in_fp) != data_size) {
                free(data);
                fprintf(stderr, "Warning: Failed to read data in block %d\n", nblocks);
                break;
            }
            
            // Compute CRC to verify integrity
            uint32_t computed_crc = crc32(0L, Z_NULL, 0);
            computed_crc = crc32(computed_crc, data, data_size);
            
            if (stored_crc != computed_crc) {
                integrity_errors++;
            }
            
            // Parse block metadata to calculate uncompressed size
            // Format: name, length, sequence, quality sections
            unsigned char *ptr = data;
            unsigned char *end = data + data_size;
            
            // Name section: u_len (4), strategy (1), c_len (4), data
            if (ptr + 9 <= end) {
                uint32_t name_usize, name_csize;
                memcpy(&name_usize, ptr, 4); ptr += 4;
                ptr++; // skip name_strat
                memcpy(&name_csize, ptr, 4); ptr += 4;
                total_uncompressed += name_usize;
                if (ptr + name_csize <= end) {
                    ptr += name_csize;
                } else {
                    // Compressed data extends beyond block - file may be truncated
                    ptr = end;
                }
            }
            
            // Length section: strategy (1), followed by variable data
            if (ptr < end) {
                uint8_t len_strat = *ptr;
                ptr++; // skip len_strat
                if (len_strat > 0) {
                    // Fixed length - skip the varint encoded length
                    int nb = 0;
                    // Simple varint decode to advance ptr
                    while (ptr + nb < end && nb < 5) {
                        if ((ptr[nb] & 0x80) == 0) {
                            nb++;
                            break;
                        }
                        nb++;
                    }
                    ptr += nb;
                } else {
                    // Variable length - skip compressed length field + data
                    if (ptr + 4 <= end) {
                        uint32_t blen;
                        memcpy(&blen, ptr, 4);
                        // Validate blen to prevent integer overflow
                        // Reasonable upper bound: blocks are typically < 1GB
                        if (blen > 0 && blen < 1000000000) {
                            // Bounds check before advancing pointer
                            if (ptr + 4 + blen <= end) {
                                ptr += 4 + blen;
                            } else {
                                ptr = end;  // Truncated data
                            }
                        } else if (blen == 0) {
                            // Empty variable length section, just skip the length field
                            ptr += 4;
                        } else {
                            // Suspicious blen value, treat as truncated
                            ptr = end;
                        }
                    }
                }
            }
            
            // Sequence section: strategy (1), u_len (4), c_len (4), data
            if (ptr + 9 <= end) {
                ptr++; // skip seq_strat
                uint32_t seq_usize, seq_csize;
                memcpy(&seq_usize, ptr, 4); ptr += 4;
                memcpy(&seq_csize, ptr, 4); ptr += 4;
                total_uncompressed += seq_usize; // Sequence bases
                if (ptr + seq_csize <= end) {
                    ptr += seq_csize;
                    
                    // Quality section: strategy (1), u_len (4), c_len (4), data
                    if (ptr + 9 <= end) {
                        ptr++; // skip qual_strat
                        uint32_t qual_usize, qual_csize;
                        memcpy(&qual_usize, ptr, 4); ptr += 4;
                        memcpy(&qual_csize, ptr, 4); ptr += 4;
                        total_uncompressed += qual_usize; // Quality scores
                    }
                } else {
                    // Compressed data extends beyond block - file may be truncated
                    ptr = end;
                }
            }
            
            // Add overhead for FASTQ format:
            // Each record has: @ + name + \n + seq + \n + + + \n + qual + \n
            // Stored names include null terminators, so overhead is:
            // 1(@) + 1(\n after seq) + 1(+) + 1(\n after +) + 1(\n after qual) = 5 bytes
            // (The \n after name replaces the stored \0)
            total_uncompressed += num_records * 5;
            
            free(data);
        } else {
            // No CRC - skip the rest of the block
            uint32_t data_size = block_size - 4;  // Subtract num_records (4)
            fseek(in_fp, data_size, SEEK_CUR);
        }
        
        nblocks++;
    }
    
    // Try to read index if present
    fqz5_index *idx = NULL;
    if (index_offset > 0) {
        idx = read_index(in_fp, index_offset);
    }
    
    // Display block information
    printf("Number of Blocks:    %d\n", nblocks);
    if (total_records > 0) {
        printf("Total Records:       %llu\n", (unsigned long long)total_records);
    }
    
    // Display size information
    if (total_uncompressed > 0) {
        printf("Uncompressed Size:   %llu bytes (%.2f MB)\n", 
               (unsigned long long)total_uncompressed, total_uncompressed / 1048576.0);
        double ratio = (double)total_uncompressed / file_size;
        printf("Compression Ratio:   %.2fx (%.2f%%)\n", ratio, (file_size * 100.0) / total_uncompressed);
    }
    
    // Try to detect interleaving
    // Note: This is a simple heuristic based on record count
    // True detection would require parsing read name patterns
    if (total_records > 0) {
        if (total_records % 2 == 0) {
            printf("Interleaved:         Possibly (even record count - heuristic)\n");
        } else {
            printf("Interleaved:         No (odd record count)\n");
        }
    }
    
    // Display index information
    if (idx) {
        printf("Index Present:       Yes (%u blocks indexed)\n", idx->nblocks);
        free_index(idx);
    } else {
        printf("Index Present:       No\n");
    }
    
    // Display integrity status
    printf("\nIntegrity Check:\n");
    if (has_crc) {
        if (integrity_errors == 0) {
            printf("  Status:            OK (all %d blocks verified)\n", nblocks);
        } else {
            printf("  Status:            FAILED (%d/%d blocks have CRC errors)\n", 
                   integrity_errors, nblocks);
        }
    } else {
        printf("  Status:            Not Available (file has no CRC checksums)\n");
        printf("  Note:              Upgrade to v1.1 format for integrity checking\n");
    }
    
    return integrity_errors > 0 ? -1 : 0;
}

// Fast integrity check - verify CRC checksums without decompressing
int check_integrity(FILE *in_fp, opts *arg) {
    uint64_t index_offset;
    int header_result = read_header(in_fp, &index_offset);
    if (header_result < 0)
        return -1;
    
    int file_version = header_result;
    
    if (file_version != 0) {
        fprintf(stderr, "Warning: File is version 1.0 or older (no CRC checksums)\n");
        fprintf(stderr, "Cannot verify integrity - file has no checksums.\n");
        return -1;
    }
    
    if (arg->verbose >= 0) {
        printf("Checking file integrity...\n");
    }
    
    int nblocks = 0;
    int errors = 0;
    
    for (;;) {
        // Check if we've reached the index
        uint64_t current_pos = ftell(in_fp);
        if (index_offset > 0 && current_pos >= index_offset)
            break;
            
        // Read block size (first 4 bytes of block)
        uint32_t block_size;
        if (fread(&block_size, 1, 4, in_fp) != 4)
            break;
        
        // Read num_records
        uint32_t num_records;
        if (fread(&num_records, 1, 4, in_fp) != 4) {
            fprintf(stderr, "ERROR: Failed to read num_records in block %d\n", nblocks);
            return -1;
        }
        
        // Read stored CRC
        uint32_t stored_crc;
        if (fread(&stored_crc, 1, 4, in_fp) != 4) {
            fprintf(stderr, "ERROR: Failed to read CRC in block %d\n", nblocks);
            return -1;
        }
        
        // Read the rest of the block data
        uint32_t data_size = block_size - 8;  // Subtract num_records (4) and CRC (4)
        unsigned char *data = malloc(data_size);
        if (!data) {
            fprintf(stderr, "ERROR: Failed to allocate memory for block %d\n", nblocks);
            return -1;
        }
        
        if (fread(data, 1, data_size, in_fp) != data_size) {
            free(data);
            fprintf(stderr, "ERROR: Failed to read data in block %d\n", nblocks);
            return -1;
        }
        
        // Compute CRC of the data
        uint32_t computed_crc = crc32(0L, Z_NULL, 0);
        computed_crc = crc32(computed_crc, data, data_size);
        free(data);
        
        nblocks++;
        
        // Verify CRC
        if (stored_crc != computed_crc) {
            fprintf(stderr, "ERROR: CRC mismatch in block %d!\n", nblocks);
            fprintf(stderr, "  Expected: 0x%08x, Got: 0x%08x\n", stored_crc, computed_crc);
            errors++;
        } else if (arg->verbose > 0) {
            printf("Block %d: CRC OK (0x%08x)\n", nblocks, stored_crc);
        }
    }
    
    if (arg->verbose >= 0) {
        if (errors == 0) {
            printf("SUCCESS: All %d blocks verified OK\n", nblocks);
        } else {
            printf("FAILED: %d/%d blocks had CRC errors\n", errors, nblocks);
        }
    }
    
    return errors > 0 ? -1 : 0;
}

void usage(FILE *fp) {
    fprintf(fp, "Usage: fqzcomp5 [options]    [input.fastq [output.fqz5]]\n");
    fprintf(fp, "Usage: fqzcomp5 [options]    [input_R1.fastq input_R2.fastq output.fqz5]\n");
    fprintf(fp, "Usage: fqzcomp5 [options] -d [input.fqz5  [output.fastq]]\n");
    fprintf(fp, "Usage: fqzcomp5 [options] -d [input.fqz5  [output_R1.fastq output_R2.fastq]]\n");
    fprintf(fp, "Usage: fqzcomp5 --check      [input.fqz5]\n");
    fprintf(fp, "Usage: fqzcomp5 --inspect    [input.fqz5]\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "    -d            Decompress\n");
    fprintf(fp, "    --check       Verify file integrity (CRC checksums) without decompressing\n");
    fprintf(fp, "    --inspect     Display comprehensive file information (version, size, ratio, etc.)\n");
    fprintf(fp, "    -p            Output name on third line (+name instead of +)\n");
    fprintf(fp, "    -t INT        Number of threads.  Defaults to 4\n");
    fprintf(fp, "    -b SIZE       Specify block size. May use K, M and G sufixes\n");
    fprintf(fp, "    -v            Increase verbostity\n");
    fprintf(fp, "    -V            Silent mode\n");
    fprintf(fp, "\n");
    fprintf(fp, "    -n INT        Name encoding method (0=rANS, 1=tok3, 2=tok3+LZP)\n");
    fprintf(fp, "    -N INT        Name encoding strategy.\n");
    fprintf(fp, "    -s INT        Sequence encoding method (0=rANS, 1=fqz)\n");
    fprintf(fp, "    -S INT        Sequence encoding strategy (context size)\n");
    fprintf(fp, "    -B            Update sequence context on both strands\n");
    fprintf(fp, "    -q INT        Quality encoding method (0=rANS, 1=fqz)\n");
    fprintf(fp, "    -Q INT        Quality encoding strategy (0 to 3)\n");
    fprintf(fp, "\n");
    fprintf(fp, "Compression levels:\n");
    fprintf(fp, "    -1            Light compression; 10MB block and rANS only\n");
    fprintf(fp, "    -3            100MB block and rANS/TOK3\n");
    fprintf(fp, "    -5            100MB block and basic seq / qual FQZ modes (default)\n");
    fprintf(fp, "    -7            500MB block and higher level FQZ modes\n");
    fprintf(fp, "    -9            Maximum compression, with 1GB blocks\n");
    fprintf(fp, "\n");
    fprintf(fp, "Paired-end mode:\n");
    fprintf(fp, "  When two input files are provided, they are automatically interleaved\n");
    fprintf(fp, "  during compression for better read-name compression.\n");
    fprintf(fp, "  When two output files are provided during decompression, the data is\n");
    fprintf(fp, "  automatically deinterleaved to separate R1 and R2 files.\n");
    fprintf(fp, "\n");
    fprintf(fp, "File Integrity:\n");
    fprintf(fp, "  FQZ5 v1.1 files include CRC32 checksums for data integrity.\n");
    fprintf(fp, "  Checksums are verified automatically during decompression.\n");
    fprintf(fp, "  Use --check to quickly verify file integrity without decompressing.\n");
    fprintf(fp, "  Use --inspect to get comprehensive information about a compressed file.\n");
}

int main(int argc, char **argv) {
    int decomp = 0;
    fqz_gparams *gp = NULL, gp_local;
    FILE *in_fp, *out_fp;
    timings t = {0};

    opts arg = {
	.qstrat = 1, // 0=rans, 1=fqz
	.qlevel = 0,
	.sstrat = 1, // 0=rans, 1=fqz
	.slevel = 12,// seq context = 4^12 
	.scustom= 0,
	.nstrat = 2, // (0=rans), 1=tok3, 2=tok3 + comments
	.nlevel = 5,
	.qauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
	         |(1<<FQZ0) |(1<<FQZ1),
	.sauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
	         |(1<<SEQ10)|(1<<SEQ12B),
	.nauto  = (1<<TLZP3)|(1<<TOK3_5_LZP),
	.both_strands =0, // adjusts seq strat 1.
	.verbose = 0,
	.blk_size = BLK_SIZE,
	.nthread = 4,
	.plus_name = 0,  // don't output name on third line by default
	.check_only = 0,  // don't do check-only mode by default
	.inspect_only = 0,  // don't do inspect-only mode by default
	.verify_crc = 1,  // verify CRC by default when available
	.paired_mode = 0,  // will be set to 1 if processing paired-end files
    };

#ifdef _WIN32
        _setmode(_fileno(stdin),  _O_BINARY);
        _setmode(_fileno(stdout), _O_BINARY);
#endif

    // Check for --check and --inspect flags before getopt processing
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--check") == 0) {
            arg.check_only = 1;
            // Remove --check from argv so it doesn't interfere with getopt
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;
            break;
        } else if (strcmp(argv[i], "--inspect") == 0) {
            arg.inspect_only = 1;
            // Remove --inspect from argv so it doesn't interfere with getopt
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;
            break;
        }
    }

    extern char *optarg;
    extern int optind;
    int opt;

    while ((opt = getopt(argc, argv, "dq:Q:b:x:Bs:S:vn:N:Vt:ph13579")) != -1) {
	switch (opt) {
	case 't':
	    arg.nthread = atoi(optarg);
	    if (arg.nthread < 1)
		arg.nthread = 1;
	    break;

	case 'v':
	    arg.verbose++;
	    break;

	case 'V':
	    arg.verbose = -1;
	    break;

	case 'd':
	    decomp = 1;
	    break;

	case 'p':
	    arg.plus_name = 1;
	    break;

	case 'B':
	    arg.both_strands=1;
	    break;

	case 's':
	    arg.sstrat = atoi(optarg);
	    if (!arg.sstrat)
		arg.sauto = 0;
	    break;
	case 'S':
	    arg.slevel = atoi(optarg);
	    arg.sstrat = 1; // -S implies -s1
	    arg.scustom = 1;
	    if (arg.slevel < 0)
		arg.slevel = 0;
	    if (arg.slevel > 16)
		arg.slevel = 16;
	    break;

	case 'n':
	    arg.nstrat = atoi(optarg);
	    arg.nauto  = 0; // override combinatorial search
	    break;
	case 'N':
	    arg.nlevel = atoi(optarg);
	    if (arg.nlevel < 0)
		arg.nlevel = 0;
	    if (arg.nlevel > 19)
		arg.nlevel = 19;
	    break;

	case 'q':
	    arg.qstrat = atoi(optarg);
	    if (arg.qstrat && !arg.qauto)
		arg.qauto = 1<<FQZ0;
	    else if (!arg.qstrat)
		arg.qauto = 0;
	    break;
	case 'Q':
	    arg.qlevel = atoi(optarg);
	    arg.qstrat = 1; // -Q implies -q1
	    arg.qauto  = 1<<(FQZ0+arg.qlevel); // override combinatorial search
	    break;

	case 'b': {
	    char *endp;
	    arg.blk_size = strtol(optarg, &endp, 0);
	    if (*endp == 'k' || *endp == 'K')
		arg.blk_size *= 1000;
	    else if (*endp == 'm' || *endp == 'M')
		arg.blk_size *= 1000000;
	    else if (*endp == 'g' || *endp == 'G')
		arg.blk_size *= 1000000000;
	    if (arg.blk_size < 1000000)
		arg.blk_size = 1000000;
	    if (arg.blk_size > 2000000000)
		arg.blk_size = 2000000000;
	    break;
	}

	case '1':
	    arg.nauto  = (1<<TLZP3);
	    arg.sauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)|(1<<LZP3);
	    arg.qauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193);
	    arg.blk_size = 10e6;
	    break;

	case '3':
	    arg.nauto  = (1<<TLZP3)|(1<<TOK3_3_LZP);
	    arg.sauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
	               | (1<<LZP3);
	    arg.qauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
		       | (1<<RANSXN1);
	    arg.blk_size = 100e6;
	    break;

	case '5':
	    arg.nauto  = (1<<TLZP3)|(1<<TOK3_5_LZP);//|(1<<TOK3_5);
	    arg.sauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
		       | (1<<LZP3) | (1<<SEQ10)|(1<<SEQ12B);
	    arg.qauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
		       | (1<<RANSXN1) | (1<<FQZ1) |(1<<FQZ3);
	    arg.blk_size = 100e6;
	    break;

	case '7':
	    arg.nauto  = (1<<TLZP3)|(1<<TOK3_7_LZP)|(1<<TOK3_7);
	    arg.sauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
		       | (1<<LZP3) |(1<<RANS65)|(1<<SEQ10)|(1<<SEQ12B)|(1<<SEQ13B);
	    arg.qauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
		       |(1<<RANS65)|(1<<FQZ0) |(1<<FQZ1)   |(1<<FQZ2)
		       |(1<<FQZ3)  |(1<<FQZ4);
	    arg.blk_size = 500e6;
	    break;

	case '9':
	    arg.nauto  = (1<<TLZP3)|(1<<TOK3_9_LZP)|(1<<TOK3_9);
	    arg.sauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
		       | (1<<RANS64)|(1<<RANS65)|(1<<RANS128)|(1<<RANS129)
		       | (1<<LZP3) | (1<<SEQ10)|(1<<SEQ12)|(1<<SEQ12B)
		       |(1<<SEQ13B)|(1<<SEQ14B);
	    arg.qauto  = (1<<RANS0)|(1<<RANS1)|(1<<RANS129)|(1<<RANS193)
		       | (1<<RANS64)|(1<<RANS65)|(1<<RANS128)|(1<<RANS129)
		       | (1<<FQZ0) |(1<<FQZ1) |(1<<FQZ2)   |(1<<FQZ3)
		       | (1<<FQZ4);
	    arg.blk_size = 1e9;
	    break;

#if 0
	case 'x': {
	    // Hex digits are:
	    // qbits  qshift
	    // pbits  pshift
	    // dbits  dshift
	    // qloc   sloc
	    // ploc   dloc
	    // do_r2  do_qavg
	    //
	    // Examples: -x 0x5570000d6e14 q40+dir =  3473340
	    //           -x 0x8252120e8d04 q4      =  724989
	    uint64_t x = strtol(optarg, NULL, 0);
	    manual_strats[manual_nstrat++] = x;

	    gp = &gp_local;
	    break;
	}
#endif

	case 'h':
	    usage(stdout);
	    return 0;

	default:
	    usage(stderr);
	    return 1;
	}
    }

    if (optind == argc && isatty(0)) {
        usage(stdout);
        return 0;
    }

    // Handle --check mode
    if (arg.check_only) {
        if (argc - optind != 1) {
            fprintf(stderr, "Error: --check requires exactly one input file\n");
            usage(stderr);
            return 1;
        }
        
        const char *in_name = argv[optind];
        FILE *in_fp = fopen(in_name, "rb");
        if (!in_fp) {
            perror(in_name);
            return 1;
        }
        
        int result = check_integrity(in_fp, &arg);
        fclose(in_fp);
        return result == 0 ? 0 : 1;
    }

    // Handle --inspect mode
    if (arg.inspect_only) {
        if (argc - optind != 1) {
            fprintf(stderr, "Error: --inspect requires exactly one input file\n");
            usage(stderr);
            return 1;
        }
        
        const char *in_name = argv[optind];
        FILE *in_fp = fopen(in_name, "rb");
        if (!in_fp) {
            perror(in_name);
            return 1;
        }
        
        int result = inspect_file(in_fp, &arg);
        fclose(in_fp);
        return result == 0 ? 0 : 1;
    }

    // Determine if we have paired files
    // Compression: 0, 1, 2, or 3 args (stdin/stdout, in, in out, in1 in2 out)
    // Decompression: 0, 1, 2, or 3 args (stdin/stdout, in, in out, in out1 out2)
    int num_file_args = argc - optind;
    int paired_mode = 0;
    
    const char *in_name1 = NULL, *in_name2 = NULL;
    const char *out_name1 = NULL, *out_name2 = NULL;
    
    if (decomp) {
        // Decompression mode
        if (num_file_args == 0) {
            // stdin to stdout
            in_name1 = NULL;
            out_name1 = NULL;
        } else if (num_file_args == 1) {
            // file to stdout
            in_name1 = argv[optind];
            out_name1 = NULL;
        } else if (num_file_args == 2) {
            // file to file
            in_name1 = argv[optind];
            out_name1 = argv[optind+1];
        } else if (num_file_args == 3) {
            // file to paired files (deinterleave)
            in_name1 = argv[optind];
            out_name1 = argv[optind+1];
            out_name2 = argv[optind+2];
            paired_mode = 1;
            arg.paired_mode = 1;  // Set in opts for encode_names
        } else {
            fprintf(stderr, "Error: Too many arguments\n");
            usage(stderr);
            return 1;
        }
    } else {
        // Compression mode
        if (num_file_args == 0) {
            // stdin to stdout
            in_name1 = NULL;
            out_name1 = NULL;
        } else if (num_file_args == 1) {
            // file to stdout
            in_name1 = argv[optind];
            out_name1 = NULL;
        } else if (num_file_args == 2) {
            // file to file
            in_name1 = argv[optind];
            out_name1 = argv[optind+1];
        } else if (num_file_args == 3) {
            // paired files to file (interleave)
            in_name1 = argv[optind];
            in_name2 = argv[optind+1];
            out_name1 = argv[optind+2];
            paired_mode = 1;
            arg.paired_mode = 1;  // Set in opts for encode_names
        } else {
            fprintf(stderr, "Error: Too many arguments\n");
            usage(stderr);
            return 1;
        }
    }
    
    // Detect if input/output should be gzipped
    int in1_is_gz = 0, in2_is_gz = 0, out1_is_gz = 0, out2_is_gz = 0;
    
    if (in_name1) {
        size_t len = strlen(in_name1);
        if (len > 3 && strcmp(in_name1 + len - 3, ".gz") == 0)
            in1_is_gz = 1;
    }
    
    if (in_name2) {
        size_t len = strlen(in_name2);
        if (len > 3 && strcmp(in_name2 + len - 3, ".gz") == 0)
            in2_is_gz = 1;
    }
    
    if (out_name1) {
        size_t len = strlen(out_name1);
        if (len > 3 && strcmp(out_name1 + len - 3, ".gz") == 0)
            out1_is_gz = 1;
    }
    
    if (out_name2) {
        size_t len = strlen(out_name2);
        if (len > 3 && strcmp(out_name2 + len - 3, ".gz") == 0)
            out2_is_gz = 1;
    }

    // Open files based on compression status
    gzFile gz_in_fp1 = NULL, gz_in_fp2 = NULL;
    gzFile gz_out_fp1 = NULL, gz_out_fp2 = NULL;
    FILE *in_fp1 = NULL, *in_fp2 = NULL;
    FILE *out_fp1 = NULL, *out_fp2 = NULL;
    
    if (decomp) {
        // For decompression: input is always .fqz5 (not gzipped)
        in_fp1 = in_name1 ? fopen(in_name1, "rb") : stdin;
        if (!in_fp1) {
            perror(in_name1);
            return 1;
        }
        
        if (paired_mode) {
            // Two output files (deinterleave)
            if (out1_is_gz) {
                gz_out_fp1 = out_name1 ? gzopen(out_name1, "wb") : gzdopen(fileno(stdout), "wb");
                if (!gz_out_fp1) {
                    perror(out_name1);
                    fclose(in_fp1);
                    return 1;
                }
            } else {
                out_fp1 = out_name1 ? fopen(out_name1, "wb") : stdout;
                if (!out_fp1) {
                    perror(out_name1);
                    fclose(in_fp1);
                    return 1;
                }
            }
            
            if (out2_is_gz) {
                gz_out_fp2 = gzopen(out_name2, "wb");
                if (!gz_out_fp2) {
                    perror(out_name2);
                    fclose(in_fp1);
                    if (out1_is_gz) gzclose(gz_out_fp1); else fclose(out_fp1);
                    return 1;
                }
            } else {
                out_fp2 = fopen(out_name2, "wb");
                if (!out_fp2) {
                    perror(out_name2);
                    fclose(in_fp1);
                    if (out1_is_gz) gzclose(gz_out_fp1); else fclose(out_fp1);
                    return 1;
                }
            }
        } else {
            // Single output file
            if (out1_is_gz) {
                gz_out_fp1 = out_name1 ? gzopen(out_name1, "wb") : gzdopen(fileno(stdout), "wb");
                if (!gz_out_fp1) {
                    perror(out_name1);
                    fclose(in_fp1);
                    return 1;
                }
            } else {
                out_fp1 = out_name1 ? fopen(out_name1, "wb") : stdout;
                if (!out_fp1) {
                    perror(out_name1);
                    fclose(in_fp1);
                    return 1;
                }
            }
        }
    } else {
        // For compression
        if (paired_mode) {
            // Two input files (interleave)
            if (in1_is_gz) {
                gz_in_fp1 = in_name1 ? gzopen(in_name1, "rb") : gzdopen(fileno(stdin), "rb");
                if (!gz_in_fp1) {
                    perror(in_name1);
                    return 1;
                }
            } else {
                // Use gzopen even for non-gzipped files (kseq handles both)
                gz_in_fp1 = in_name1 ? gzopen(in_name1, "rb") : gzdopen(fileno(stdin), "rb");
                if (!gz_in_fp1) {
                    perror(in_name1);
                    return 1;
                }
            }
            
            if (in2_is_gz) {
                gz_in_fp2 = gzopen(in_name2, "rb");
                if (!gz_in_fp2) {
                    perror(in_name2);
                    if (in1_is_gz) gzclose(gz_in_fp1); else fclose(in_fp1);
                    return 1;
                }
            } else {
                // Use gzopen even for non-gzipped files (kseq handles both)
                gz_in_fp2 = gzopen(in_name2, "rb");
                if (!gz_in_fp2) {
                    perror(in_name2);
                    gzclose(gz_in_fp1);
                    return 1;
                }
            }
        } else {
            // Single input file
            if (in1_is_gz) {
                gz_in_fp1 = in_name1 ? gzopen(in_name1, "rb") : gzdopen(fileno(stdin), "rb");
                if (!gz_in_fp1) {
                    perror(in_name1);
                    return 1;
                }
            } else {
                // Use gzopen even for non-gzipped files (kseq handles both)
                gz_in_fp1 = in_name1 ? gzopen(in_name1, "rb") : gzdopen(fileno(stdin), "rb");
                if (!gz_in_fp1) {
                    perror(in_name1);
                    return 1;
                }
            }
        }
        
        out_fp1 = out_name1 ? fopen(out_name1, "wb") : stdout;
        if (!out_fp1) {
            perror(out_name1);
            if (paired_mode) {
                if (in1_is_gz) gzclose(gz_in_fp1); else fclose(in_fp1);
                if (in2_is_gz) gzclose(gz_in_fp2); else fclose(in_fp2);
            } else {
                if (in1_is_gz) gzclose(gz_in_fp1); else fclose(in_fp1);
            }
            return 1;
        }
    }

    // Block based, for arbitrary sizes of input
    if (decomp) {
        if (paired_mode) {
            // Deinterleave to two files
            if (out1_is_gz && out2_is_gz) {
                if (decode_gzip_deinterleaved(in_fp1, gz_out_fp1, gz_out_fp2, &arg, &t) < 0)
                    exit(1);
            } else if (!out1_is_gz && !out2_is_gz) {
                if (decode_deinterleaved(in_fp1, out_fp1, out_fp2, &arg, &t) < 0)
                    exit(1);
            } else {
                fprintf(stderr, "Error: Both output files must have the same format (both .gz or both plain)\n");
                exit(1);
            }
        } else {
            // Single output file
            if (out1_is_gz) {
                if (decode_gzip(in_fp1, gz_out_fp1, &arg, &t) < 0)
                    exit(1);
            } else {
                if (decode(in_fp1, out_fp1, &arg, &t) < 0)
                    exit(1);
            }
        }
    } else {
        if (paired_mode) {
            // Interleave two input files
            // Both inputs should be gzipped or both plain
            // Both files are now opened with gzopen
            if (encode_interleaved(gz_in_fp1, gz_in_fp2, out_fp1, gp, &arg, &t) < 0)
                exit(1);
        } else {
            // Single input file
            if (encode_gzip(gz_in_fp1, out_fp1, gp, &arg, &t) < 0)
                exit(1);
        }
    }

    if (arg.verbose >= 0) {
        fprintf(stderr, "All %ld blocks combined:\n", t.nblock);
        fprintf(stderr, "Names    %10ld to %10ld in %.2f sec\n",
                t.nusize, t.ncsize, t.ntime/1e6);
        fprintf(stderr, "Lengths  %10ld to %10ld\n",
                t.lusize, t.lcsize);
        fprintf(stderr, "Seqs     %10ld to %10ld in %.2f sec\n", 
                t.susize, t.scsize, t.stime/1e6);
        fprintf(stderr, "Qual     %10ld to %10ld in %.2f sec\n", 
                t.qusize, t.qcsize, t.qtime/1e6);
    }

    // Close files
    if (decomp) {
        fclose(in_fp1);
        if (paired_mode) {
            if (out1_is_gz) gzclose(gz_out_fp1); else fclose(out_fp1);
            if (out2_is_gz) gzclose(gz_out_fp2); else fclose(out_fp2);
        } else {
            if (out1_is_gz) gzclose(gz_out_fp1); else fclose(out_fp1);
        }
    } else {
        if (paired_mode) {
            gzclose(gz_in_fp1);
            gzclose(gz_in_fp2);
        } else {
            gzclose(gz_in_fp1);
        }
        fclose(out_fp1);
    }

    return 0;
}
