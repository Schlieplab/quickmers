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
    int variable_length;  /* 1: results carry a length sentinel bit */

    /* Operation strings */
    KmerSet* ops_set;
    KmerSet* ops_cur;

    /* Deduplication */
    BitKmerSet* seen;

    /* Current batch */
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

uint64_t fixed_length_levenshtein_ball_iter_next(LevBallIter* it,
                                                 int* has_value);

void fixed_length_levenshtein_ball_iter_free(LevBallIter* it);

/* ---------------- Variable-length List API ----------------
 * Results are sentinel-encoded: a 1 bit sits directly above the 2*len
 * bits of the k-mer, so k-mers of different lengths never collide
 * (e.g. "AC" -> 0b10001, "AAC" -> 0b1000001). Use
 * variable_kmer_length / variable_kmer_value / decode_variable_kmer.
 * Requires strlen(kmer) + radius <= LEV_BALL_MAX_VAR_LEN.
 */
#define LEV_BALL_MAX_VAR_LEN 31

void levenshtein_ball(const char* kmer,
                      int radius,
                      uint64_t** out,
                      size_t* out_count);

/* ---------------- Variable-length Iterator API ---------------- */
void levenshtein_ball_iter_init(LevBallIter* it,
                                const char* kmer,
                                int radius);

uint64_t levenshtein_ball_iter_next(LevBallIter* it,
                                    int* has_value);

void levenshtein_ball_iter_free(LevBallIter* it);

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

/* Length of a sentinel-encoded k-mer (position of the sentinel bit / 2) */
static inline int variable_kmer_length(uint64_t x)
{
    int len = 0;
    while ((x >> (2 * len)) > 1) len++;
    return len;
}

/* Sentinel-encoded k-mer -> plain encoding (same as encode_kmer) */
static inline uint64_t variable_kmer_value(uint64_t x, int len)
{
    return x ^ ((uint64_t)1 << (2 * len));
}

/* Decode a sentinel-encoded k-mer; out needs LEV_BALL_MAX_VAR_LEN + 1 bytes */
static inline void decode_variable_kmer(uint64_t x, char* out)
{
    int len = variable_kmer_length(x);
    decode_kmer(variable_kmer_value(x, len), len, out);
}

#endif // LEVENSHTEIN_BALL_H