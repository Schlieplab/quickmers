/*
Filename: quickmers/include/levenshtein.h
Author: Kian Jalilian
Copyright: 2025, Kian Jalilian
Version: 0.1.0
Description: Header file for levenshtein distance implementation in C
License: LGPL-3.0-or-later
*/
#ifndef LEVENSHTEIN_H
#define LEVENSHTEIN_H

#include <stdint.h>

int64_t myers(uint8_t *t, int64_t n, uint8_t *p, int64_t m);
void myers_batch_avx2(
    const uint8_t *query, int64_t qlen,
    const uint8_t **kmers, int64_t n_kmers, int64_t kmer_len,
    int64_t *out
);
void myers_dispatch(
    const uint8_t *query, int64_t qlen,
    const uint8_t **kmers, int64_t n_kmers, int64_t kmer_len,
    int64_t *out
);
#endif