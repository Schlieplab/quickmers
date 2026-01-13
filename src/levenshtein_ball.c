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

static const char ALPHABET[4] = {'A','C','G','T'};

/* ============================================================================
 *                          Hash-set Utilities
 * ============================================================================
 */

/*
 * Insert a string into a KmerSet (hash set of strings).
 *
 * Parameters:
 *   set  - Pointer to the hash set.
 *   kmer - String to insert.
 *
 * Behavior:
 *   - If the string is not already present, it is inserted.
 *   - If the string already exists, the passed-in string is freed to avoid
 *     duplicates.
 */
static void kmerset_add(KmerSet** set, char* kmer) {
    KmerSet* s = NULL;
    HASH_FIND_STR(*set, kmer, s);
    if (!s) {
        s = (KmerSet*)malloc(sizeof(KmerSet));
        s->kmer = kmer;
        HASH_ADD_KEYPTR(hh, *set, s->kmer, strlen(s->kmer), s);
    } else {
        free(kmer); // duplicate, discard
    }
}

/*
 * Convert a KmerSet into a dynamically allocated array of strings.
 *
 * Parameters:
 *   set       - Hash set containing k-mer strings.
 *   out_kmers - Output pointer to array of strings.
 *   out_count - Output number of strings in the array.
 *
 * Behavior:
 *   - Allocates a char** array and transfers ownership of the strings.
 *   - Frees the hash table nodes, but not the strings themselves.
 */
static void kmerset_to_array(KmerSet* set, char*** out_kmers, size_t* out_count) {
    size_t count = HASH_COUNT(set);
    char** arr = (char**)malloc(sizeof(char*) * count);
    size_t idx = 0;
    KmerSet* s;
    KmerSet* tmp;
    HASH_ITER(hh, set, s, tmp) {
        arr[idx++] = s->kmer;
        HASH_DEL(set, s);
        free(s); // free struct but keep string
    }
    *out_kmers = arr;
    *out_count = count;
}

/* ============================================================================
 *                    Unique Permutations of Operation Strings
 * ============================================================================
 */

/*
 * Helper structure for generating unique permutations of a multiset of
 * edit operations (e.g., M, C, I, D).
 *
 * Fields:
 *   ops    - Array of distinct operation symbols.
 *   counts - Remaining counts for each operation symbol.
 *   n      - Number of unique operation symbols.
 *   length - Total length of the permutation to generate.
 */
typedef struct {
    char ops[64];
    int counts[4];
    int n;
    int length;
} OpCounter;

/*
 * Recursively generate all unique permutations of the multiset encoded
 * in the OpCounter and insert them into a hash set.
 *
 * Parameters:
 *   ctr      - Pointer to OpCounter describing available symbols and counts.
 *   buffer   - Temporary buffer used to build the permutation.
 *   depth    - Current recursion depth.
 *   out_set  - Hash set that will store unique operation strings.
 *
 * Behavior:
 *   - When depth == ctr->length, a complete operation string is formed.
 *   - Each unique permutation is inserted into out_set.
 */
static void unique_permutations(OpCounter* ctr, char* buffer, int depth, KmerSet** out_set) {
    if (depth == ctr->length) {
        char* s = (char*)malloc(ctr->length + 1);
        memcpy(s, buffer, ctr->length);
        s[ctr->length] = '\0';
        kmerset_add(out_set, s);
        return;
    }
    for (int i = 0; i < ctr->n; i++) {
        if (ctr->counts[i] > 0) {
            ctr->counts[i]--;
            buffer[depth] = ctr->ops[i];
            unique_permutations(ctr, buffer, depth + 1, out_set);
            ctr->counts[i]++;
        }
    }
}

/* ============================================================================
 *                   Generation of Valid Operation Strings
 * ============================================================================
 */

/*
 * Generate all valid edit operation strings for a fixed-length Levenshtein ball.
 *
 * Parameters:
 *   k        - Length of the input k-mer.
 *   radius   - Maximum allowed edit distance.
 *   out_set  - Hash set to store generated operation strings.
 *
 * Operation symbols:
 *   M = Match (copy character)
 *   C = Change (substitution)
 *   I = Insert (add new character)
 *   D = Delete (skip character)
 *
 * Constraints enforced:
 *   - The total edit cost must be <= radius.
 *   - The final string length must equal the original length.
 *     This is enforced by requiring #I == #D.
 *
 * Logic:
 *   - x = number of substitutions (C)
 *   - y = number of insertion/deletion pairs (I and D)
 *   - Number of matches M is computed so total operations consume k input chars.
 *   - All unique permutations of these operations are generated.
 */
static void generate_operation_strings(int k, int radius, KmerSet** out_set) {
    for (int x = 0; x <= radius; x++) {
        for (int y = 0; y <= (radius - x)/2; y++) {
            int num_C = x, num_I = y, num_D = y, num_M = k - x - y;
            if (num_M < 0) continue;

            /* Build a multiset of operations */
            char ops[64];
            int idx = 0;
            for (int i = 0; i < num_M; i++) ops[idx++] = 'M';
            for (int i = 0; i < num_C; i++) ops[idx++] = 'C';
            for (int i = 0; i < num_I; i++) ops[idx++] = 'I';
            for (int i = 0; i < num_D; i++) ops[idx++] = 'D';

            /* Count unique symbols */
            char unique_ops[4];
            int counts[4] = {0};
            int n_unique = 0;
            for (int i = 0; i < idx; i++) {
                char c = ops[i];
                int found = -1;
                for (int j = 0; j < n_unique; j++) {
                    if (unique_ops[j] == c) { found = j; break; }
                }
                if (found >= 0) counts[found]++;
                else { unique_ops[n_unique] = c; counts[n_unique] = 1; n_unique++; }
            }

            /* Generate all unique permutations */
            OpCounter ctr = {.n = n_unique, .length = idx};
            memcpy(ctr.ops, unique_ops, n_unique);
            memcpy(ctr.counts, counts, n_unique * sizeof(int));

            char buffer[64];
            unique_permutations(&ctr, buffer, 0, out_set);
        }
    }
}

/* ============================================================================
 *            Applying Operations (Full Materialization Version)
 * ============================================================================
 */

/*
 * Apply a given operation string to the input k-mer and insert all resulting
 * strings into a hash set.
 *
 * Parameters:
 *   kmer     - Input string.
 *   opstring - Sequence of operations (M, C, I, D).
 *   out_set  - Hash set collecting all generated results.
 *
 * Behavior:
 *   - Maintains a dynamically growing list of partial results.
 *   - Each operation transforms the current set of strings.
 *   - At the end, all results are inserted into out_set with deduplication.
 */
static void apply_operations(const char* kmer, const char* opstring, KmerSet** out_set) {
    char** results = (char**)malloc(sizeof(char*));
    results[0] = strdup("");
    size_t res_count = 1;

    int i = 0;  /* index into original k-mer */
    for (const char* op = opstring; *op; op++) {
        size_t new_count = 0;
        size_t cap = res_count * 4;
        char** new_results = (char**)malloc(sizeof(char*) * cap);
        for (size_t r = 0; r < res_count; r++) {
            char* prev = results[r];
            size_t len = strlen(prev);
            if (*op == 'M') {
                char* s = (char*)malloc(len + 2);
                memcpy(s, prev, len); s[len] = kmer[i]; s[len+1] = '\0';
                new_results[new_count++] = s;
            } else if (*op == 'C') {
                for (int b = 0; b < 4; b++) {
                    if (ALPHABET[b] != kmer[i]) {
                        char* s = (char*)malloc(len + 2);
                        memcpy(s, prev, len); s[len] = ALPHABET[b]; s[len+1] = '\0';
                        new_results[new_count++] = s;
                    }
                }
            } else if (*op == 'I') {
                for (int b = 0; b < 4; b++) {
                    char* s = (char*)malloc(len + 2);
                    memcpy(s, prev, len); s[len] = ALPHABET[b]; s[len+1] = '\0';
                    new_results[new_count++] = s;
                }
            } else if (*op == 'D') {
                char* s = strdup(prev);
                new_results[new_count++] = s;
                continue;
            }
            free(prev);
        }
        free(results);
        results = new_results;
        res_count = new_count;
        i += (*op == 'M' || *op == 'C' || *op == 'D');
    }

    /* Insert all generated strings into the output set */
    for (size_t r = 0; r < res_count; r++) {
        kmerset_add(out_set, results[r]);
    }
    free(results);
}

/* ============================================================================
 *          Applying Operations (Collector for Iterator API)
 * ============================================================================
 */

/*
 * Apply a given operation string to the input k-mer and collect all resulting
 * strings into an array (without deduplication).
 *
 * Parameters:
 *   kmer        - Input string.
 *   opstring    - Sequence of operations.
 *   out_results - Output pointer to array of generated strings.
 *   out_count   - Output number of generated strings.
 *
 * Behavior:
 *   - Similar to apply_operations(), but instead of inserting into a hash set,
 *     it returns the full batch for iteration.
 *   - Caller is responsible for freeing the returned array and strings.
 */
static void apply_operations_collect(const char* kmer, const char* opstring,
                                     char*** out_results, size_t* out_count)
{
    char** results = (char**)malloc(sizeof(char*));
    results[0] = strdup("");
    size_t res_count = 1;

    int i = 0;
    for (const char* op = opstring; *op; op++) {
        size_t new_count = 0;
        size_t cap = res_count * 4;
        char** new_results = (char**)malloc(sizeof(char*) * cap);
        for (size_t r = 0; r < res_count; r++) {
            char* prev = results[r];
            size_t len = strlen(prev);
            if (*op == 'M') {
                char* s = (char*)malloc(len + 2);
                memcpy(s, prev, len); s[len] = kmer[i]; s[len+1] = '\0';
                new_results[new_count++] = s;
            } else if (*op == 'C') {
                for (int b = 0; b < 4; b++) {
                    if (ALPHABET[b] != kmer[i]) {
                        char* s = (char*)malloc(len + 2);
                        memcpy(s, prev, len); s[len] = ALPHABET[b]; s[len+1] = '\0';
                        new_results[new_count++] = s;
                    }
                }
            } else if (*op == 'I') {
                for (int b = 0; b < 4; b++) {
                    char* s = (char*)malloc(len + 2);
                    memcpy(s, prev, len); s[len] = ALPHABET[b]; s[len+1] = '\0';
                    new_results[new_count++] = s;
                }
            } else if (*op == 'D') {
                char* s = strdup(prev);
                new_results[new_count++] = s;
                continue;
            }
            free(prev);
        }
        free(results);
        results = new_results;
        res_count = new_count;
        i += (*op == 'M' || *op == 'C' || *op == 'D');
    }

    *out_results = results;
    *out_count = res_count;
}

/* ============================================================================
 *                             Iterator API
 * ============================================================================
 */

/*
 * Initialize a Levenshtein ball iterator.
 *
 * Parameters:
 *   it     - Pointer to iterator structure.
 *   kmer   - Input string.
 *   radius - Maximum allowed edit distance.
 *
 * Behavior:
 *   - Precomputes all valid operation strings.
 *   - Initializes internal state for lazy iteration.
 */
void fixed_length_levenshtein_ball_iter_init(LevBallIter* it,
                                            const char* kmer,
                                            int radius)
{
    memset(it, 0, sizeof(*it));
    it->kmer = strdup(kmer);
    it->radius = radius;

    generate_operation_strings((int)strlen(kmer), radius, &it->ops_set);
    it->ops_cur = it->ops_set;
    it->seen = NULL;
    it->results = NULL;
    it->res_count = 0;
    it->res_index = 0;
    it->initialized = 1;
}

/*
 * Return the next unique k-mer in the fixed-length Levenshtein ball.
 *
 * Parameters:
 *   it - Pointer to an initialized iterator.
 *
 * Returns:
 *   - A newly allocated string (caller takes ownership), or
 *   - NULL when all results have been exhausted.
 *
 * Behavior:
 *   - Iterates over batches of results generated from each operation string.
 *   - Deduplicates results across all batches using a hash set.
 */
char* fixed_length_levenshtein_ball_iter_next(LevBallIter* it)
{
    if (!it->initialized) return NULL;

    while (1) {

        /* 1. Return next result from current batch, if available */
        if (it->results && it->res_index < it->res_count) {
            char* s = it->results[it->res_index++];

            /* Deduplicate across all previously returned strings */
            KmerSet* found = NULL;
            HASH_FIND_STR(it->seen, s, found);
            if (!found) {
                kmerset_add(&it->seen, strdup(s));
                return s;   // caller takes ownership
            } else {
                free(s);
                continue;
            }
        }

        /* 2. Free the previous batch */
        if (it->results) {
            free(it->results);
            it->results = NULL;
        }

        /* 3. Move to the next operation string */
        if (!it->ops_cur) return NULL;  // fully exhausted

        char* op = it->ops_cur->kmer;
        it->ops_cur = it->ops_cur->hh.next;

        /* 4. Generate the next batch of results */
        apply_operations_collect(it->kmer, op,
                                 &it->results, &it->res_count);
        it->res_index = 0;
    }
}

/*
 * Free all resources associated with a Levenshtein ball iterator.
 *
 * Parameters:
 *   it - Pointer to iterator structure.
 *
 * Behavior:
 *   - Frees all internally allocated strings, hash tables, and buffers.
 */
void fixed_length_levenshtein_ball_iter_free(LevBallIter* it)
{
    KmerSet* s;
    KmerSet* tmp;

    /* Free all seen result strings */
    HASH_ITER(hh, it->seen, s, tmp) {
        free(s->kmer);
        HASH_DEL(it->seen, s);
        free(s);
    }

    /* Free all generated operation strings */
    HASH_ITER(hh, it->ops_set, s, tmp) {
        free(s->kmer);
        HASH_DEL(it->ops_set, s);
        free(s);
    }

    /* Free any remaining batch results */
    if (it->results) {
        for (size_t i = 0; i < it->res_count; i++) free(it->results[i]);
        free(it->results);
    }

    free((char*)it->kmer);
}

/* ============================================================================
 *                             List API
 * ============================================================================
 */

/*
 * Generate the full fixed-length Levenshtein ball and return it as a list.
 *
 * Parameters:
 *   kmer       - Input string.
 *   radius     - Maximum allowed edit distance.
 *   out_kmers  - Output pointer to array of result strings.
 *   out_count  - Output number of result strings.
 *
 * Behavior:
 *   - Eagerly generates all valid neighbors.
 *   - Deduplicates results using a hash set.
 *   - Returns a dynamically allocated array of strings.
 *
 * Ownership:
 *   - Caller owns the returned array and all strings within it.
 */
void fixed_length_levenshtein_ball(const char* kmer, int radius,
                                   char*** out_kmers, size_t* out_count)
{
    KmerSet* ops_set = NULL;
    generate_operation_strings((int)strlen(kmer), radius, &ops_set);

    KmerSet* all_set = NULL;
    KmerSet* node;
    KmerSet* tmp;

    /* Apply each operation string and collect all results */
    HASH_ITER(hh, ops_set, node, tmp) {
        apply_operations(kmer, node->kmer, &all_set);
        free(node->kmer);
        HASH_DEL(ops_set, node);
        free(node);
    }

    /* Convert final set to array */
    kmerset_to_array(all_set, out_kmers, out_count);
}
