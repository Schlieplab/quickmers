#ifndef LEVENSHTEIN_BALL_H
#define LEVENSHTEIN_BALL_H

#include <stddef.h>
#include "uthash.h"

/* Forward declare KmerSet if it's defined elsewhere */
typedef struct {
    char* kmer;
    UT_hash_handle hh;
} KmerSet;

/* ---------------- Used for the iterator ---------------- */
typedef struct {
    const char* kmer;
    int radius;

    /* Operation strings */
    KmerSet* ops_set;
    KmerSet* ops_cur;

    /* Current apply_operations state */
    char** results;
    size_t res_count;
    size_t res_index;

    /* Deduplication */
    KmerSet* seen;

    /* Flags */
    int initialized;
} LevBallIter;

/* List API */
void fixed_length_levenshtein_ball(const char* kmer, int radius,
                                   char*** out_kmers, size_t* out_count);

/* Iterator API */
void fixed_length_levenshtein_ball_iter_init(LevBallIter* it,
                                             const char* kmer,
                                             int radius);

char* fixed_length_levenshtein_ball_iter_next(LevBallIter* it);

void fixed_length_levenshtein_ball_iter_free(LevBallIter* it);

#endif // LEVENSHTEIN_BALL_H
