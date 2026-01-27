#ifndef LEVENSHTEIN_BALL_H
#define LEVENSHTEIN_BALL_H

#include <stddef.h>
#include <stdint.h>
#include "uthash.h"

typedef struct {
    char* kmer;
    UT_hash_handle hh;
} KmerSet;

/* ---------------- Hash set for uint64_t kmers ---------------- */
typedef struct BitKmerSet {
    uint64_t key;
    UT_hash_handle hh;
} BitKmerSet;

/* ---------------- Iterator struct ---------------- */
typedef struct {
    int k;                /* length of k-mer */
    uint64_t encoded;     /* input k-mer encoded as uint64_t */

    /* Operation strings (from MCID generation) */
    KmerSet* ops_set;     /* set of opstrings */
    KmerSet* ops_cur;     /* current opstring being iterated */

    /* Deduplication of returned values */
    BitKmerSet* seen;

    /* Current batch state (optional, can be NULL) */
    uint64_t* results;
    size_t res_count;
    size_t res_index;
} LevBallIter;

/* ---------------- List API ---------------- */
void fixed_length_levenshtein_ball(const char* kmer,
                                   int radius,
                                   uint64_t** out,
                                   size_t* out_count);

/* ---------------- Iterator API ---------------- */
void fixed_length_levenshtein_ball_iter_init(LevBallIter* it,
                                             const char* kmer,
                                             int radius);

/* Returns next value; sets *has_value=0 when exhausted */
uint64_t fixed_length_levenshtein_ball_iter_next(LevBallIter* it,
                                                 int* has_value);

void fixed_length_levenshtein_ball_iter_free(LevBallIter* it);

/* ---------------- Encoding helpers ---------------- */
static inline uint64_t encode_kmer(const char* kmer, int k)
{
    uint64_t x = 0;
    for (int i = 0; i < k; i++) {
        x <<= 2;
        switch (kmer[i]) {
            case 'A': x |= 0; break;
            case 'C': x |= 1; break;
            case 'G': x |= 2; break;
            case 'T': x |= 3; break;
        }
    }
    return x;
}

static inline void decode_kmer(uint64_t x, int k, char* out)
{
    for (int i = k - 1; i >= 0; i--) {
        out[i] = "ACGT"[x & 3];
        x >>= 2;
    }
    out[k] = '\0';
}

#endif // LEVENSHTEIN_BALL_H