/*
Filename: src/levenshtein_ball.c
Author: Kian Jalilian
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: Fixed-length Levenshtein ball generation using uthash
License: LGPL-3.0-or-later

This file implements generation of the fixed-length Levenshtein ball for
DNA k-mers over the alphabet {A, C, G, T}. It supports two APIs:

1. List API:
   - Generates all neighbors eagerly and returns them as a list.

2. Iterator API:
   - Lazily generates neighbors one at a time, suitable for bindings
     (e.g., Python generators).

The core idea is to:
- Enumerate all valid edit operation sequences ("M", "C", "I", "D") whose
  total edit cost is <= radius and whose net effect preserves string length.
- Apply each operation sequence to the input k-mer.
- Deduplicate results using hash sets (uthash).
*/

#include "levenshtein_ball.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ============================================================================
 *                          Hash set (uint64_t)
 * ============================================================================ */
static inline void bitset_add(BitKmerSet** set, uint64_t key)
{
    BitKmerSet* s;
    HASH_FIND(hh, *set, &key, sizeof(uint64_t), s);
    if (!s) {
        s = malloc(sizeof(*s));
        s->key = key;
        HASH_ADD(hh, *set, key, sizeof(uint64_t), s);
    }
}

/* Add a kmer string to a KmerSet hash table */
static inline void kmerset_add(KmerSet** set, const char* kmer) {
    KmerSet* entry = NULL;

    /* Check if it already exists */
    HASH_FIND_STR(*set, kmer, entry);
    if (entry) {
        return; // already present
    }

    /* Allocate new entry */
    entry = (KmerSet*)malloc(sizeof(KmerSet));
    if (!entry) {
        return; // malloc failed, could handle error better
    }

    entry->kmer = strdup(kmer); // copy string
    HASH_ADD_KEYPTR(hh, *set, entry->kmer, strlen(entry->kmer), entry);
}

/* ============================================================================
 *                   Operation String Generation (unchanged)
 * ============================================================================ */

typedef struct {
    uint8_t ops[64];
    int counts[4];
    int n;
    int length;
} OpCounter;

static void unique_permutations(OpCounter* ctr,
                                char* buffer,
                                int depth,
                                KmerSet** out)
{
    if (depth == ctr->length) {
        char* s = malloc(ctr->length + 1);
        memcpy(s, buffer, ctr->length);
        s[ctr->length] = '\0';
        kmerset_add(out, s);
        return;
    }
    for (int i = 0; i < ctr->n; i++) {
        if (ctr->counts[i] > 0) {
            ctr->counts[i]--;
            buffer[depth] = ctr->ops[i];
            unique_permutations(ctr, buffer, depth + 1, out);
            ctr->counts[i]++;
        }
    }
}

static void generate_operation_strings(int k, int radius, KmerSet** out)
{
    for (int y = 0; y <= radius / 2; y++) {
        int x = radius - 2 * y;  // total substitutions minus 2*indels

        int M = k - x - y;
        if (M < 0) continue;

        char ops[64];
        int idx = 0;
        for (int i = 0; i < M; i++) ops[idx++] = 'M';
        for (int i = 0; i < x; i++) ops[idx++] = 'C';
        for (int i = 0; i < y; i++) ops[idx++] = 'I';
        for (int i = 0; i < y; i++) ops[idx++] = 'D';

        // Build unique counts for permutations
        char uniq[4];
        int cnt[4] = {0};
        int n = 0;

        for (int i = 0; i < idx; i++) {
            int j;
            for (j = 0; j < n; j++)
                if (uniq[j] == ops[i]) break;
            if (j == n) uniq[n++] = ops[i];
            cnt[j]++;
        }

        OpCounter ctr = {.n = n, .length = idx};
        memcpy(ctr.ops, uniq, n);
        memcpy(ctr.counts, cnt, n * sizeof(int));

        char buf[64];
        unique_permutations(&ctr, buf, 0, out);
    }
}


/* ============================================================================
 *                   Apply Operations (binary)
 * ============================================================================ */

static void apply_operations_binary(uint64_t encoded,
                                    int k,
                                    const char* ops,
                                    BitKmerSet** out)
{
    uint64_t* states = malloc(sizeof(uint64_t));
    size_t n = 1;
    states[0] = 0;

    int src = 0;

    for (const char* op = ops; *op; op++) {
        uint64_t* next = malloc(sizeof(uint64_t) * n * 4);
        size_t next_n = 0;

        for (size_t i = 0; i < n; i++) {
            uint64_t cur = states[i];

            if (*op == 'M') {
                uint64_t b = (encoded >> (2 * (k - src - 1))) & 3;
                next[next_n++] = (cur << 2) | b;
            } else if (*op == 'C') {
                uint64_t b0 = (encoded >> (2 * (k - src - 1))) & 3;
                for (uint64_t b = 0; b < 4; b++)
                    if (b != b0)
                        next[next_n++] = (cur << 2) | b;
            } else if (*op == 'I') {
                for (uint64_t b = 0; b < 4; b++)
                    next[next_n++] = (cur << 2) | b;
            } else if (*op == 'D') {
                next[next_n++] = cur;
            }
        }

        free(states);
        states = next;
        n = next_n;

        if (*op != 'I') src++;
    }

    for (size_t i = 0; i < n; i++)
        bitset_add(out, states[i]);

    free(states);
}

/* ============================================================================
 *                             List API
 * ============================================================================ */

void fixed_length_levenshtein_ball(const char* kmer,
                                   int radius,
                                   uint64_t** out,
                                   size_t* out_count)
{
    int k = strlen(kmer);
    uint64_t encoded = encode_kmer(kmer, k);

    KmerSet* ops = NULL;
    generate_operation_strings(k, radius, &ops);

    BitKmerSet* result = NULL;

    KmerSet* s;
    KmerSet* tmp;
    HASH_ITER(hh, ops, s, tmp) {
        apply_operations_binary(encoded, k, s->kmer, &result);
        free(s->kmer);
        HASH_DEL(ops, s);
        free(s);
    }

    size_t n = HASH_COUNT(result);
    uint64_t* arr = malloc(sizeof(uint64_t) * n);

    size_t i = 0;
    BitKmerSet* b;
    BitKmerSet* btmp;
    HASH_ITER(hh, result, b, btmp) {
        arr[i++] = b->key;
        HASH_DEL(result, b);
        free(b);
    }

    *out = arr;
    *out_count = n;
}

/* ============================================================================
 *                             Iterator API
 * ============================================================================ */

void fixed_length_levenshtein_ball_iter_init(LevBallIter* it,
                                             const char* kmer,
                                             int radius)
{
    memset(it, 0, sizeof(*it));
    it->k = strlen(kmer);
    it->encoded = encode_kmer(kmer, it->k);

    generate_operation_strings(it->k, radius, &it->ops_set);
    it->ops_cur = it->ops_set;
}

uint64_t fixed_length_levenshtein_ball_iter_next(LevBallIter* it,
                                                 int* has_value)
{
    while (1) {
        // 1. Return from current results array if not exhausted
        if (it->results && it->res_index < it->res_count) {
            uint64_t v = it->results[it->res_index++];

            // Deduplicate
            BitKmerSet* f;
            HASH_FIND(hh, it->seen, &v, sizeof(uint64_t), f);
            if (!f) {
                bitset_add(&it->seen, v);
                *has_value = 1;
                return v;
            }
            // Already seen, continue to next in results
            continue;
        }

        // 2. Free previous results array
        if (it->results) {
            free(it->results);
            it->results = NULL;
            it->res_count = 0;
            it->res_index = 0;
        }

        // 3. Check if any operation strings remain
        if (!it->ops_cur) {
            *has_value = 0;
            return 0;
        }

        // 4. Take next operation string
        char* op = it->ops_cur->kmer;
        it->ops_cur = it->ops_cur->hh.next;

        // 5. Apply operations to generate a temporary set
        BitKmerSet* tmp_set = NULL;
        apply_operations_binary(it->encoded, it->k, op, &tmp_set);

        // 6. Convert tmp_set to array for yielding
        size_t n = HASH_COUNT(tmp_set);
        if (n == 0) continue; // skip empty results

        it->results = malloc(sizeof(uint64_t) * n);
        it->res_count = n;
        it->res_index = 0;

        size_t i = 0;
        BitKmerSet* b;
        BitKmerSet* bt;
        HASH_ITER(hh, tmp_set, b, bt) {
            it->results[i++] = b->key;
            HASH_DEL(tmp_set, b);
            free(b);
        }

        // Loop will now go back to the top and return first value from results
    }
}


void fixed_length_levenshtein_ball_iter_free(LevBallIter* it)
{
    BitKmerSet* b;
    BitKmerSet* btmp;
    HASH_ITER(hh, it->seen, b, btmp) {
        HASH_DEL(it->seen, b);
        free(b);
    }

    KmerSet* s;
    KmerSet* stmp;
    HASH_ITER(hh, it->ops_set, s, stmp) {
        free(s->kmer);
        HASH_DEL(it->ops_set, s);
        free(s);
    }
}
