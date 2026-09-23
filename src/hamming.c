/*
Filename: src/hamming.c
Author: Kian Jalilian
Copyright: 2025, Alexander Schliep
Version: 0.2.0
Description: Function implementations for hamming distance calculations
License: LGPL-3.0-or-later
*/
#include "hamming.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static const uint8_t base_to_bits[256] = {
    ['A'] = 0b00,
    ['C'] = 0b01,
    ['G'] = 0b10,
    ['T'] = 0b11,
};

uint32_t kmer_to_uint32(const char *s, int k) {
    uint32_t result = 0;
    for (int i = 0; i < k; i++) {
        uint8_t bits = base_to_bits[(unsigned char)s[i]];
        result = (result << 2) | bits;
    }
    return result;
}

uint64_t kmer_to_uint64(const char *s, int k) {
    uint64_t result = 0;
    for (int i = 0; i < k; i++) {
        uint8_t bits = base_to_bits[(unsigned char)s[i]];
        result = (result << 2) | bits;
    }
    return result;
}

void hamming_distance_array_32bit(
    const char *query,
    const char **kmers_list,
    size_t n,
    int *distances)
{
    int k = strlen(query); 

    // Encode query k-mer once
    uint32_t query_kmer = kmer_to_uint32(query, k);

    // Allocate temporary array for all kmers encoded
    uint32_t* kmers_encoded = (uint32_t*)malloc(n * sizeof(uint32_t));
    if (!kmers_encoded) return;  // allocation failure

    // Encode all kmers first
    for (size_t i = 0; i < n; i++) {
        // Make sure the input string is valid
        if (!kmers_list[i]) {
            kmers_encoded[i] = 0;  // or handle differently
        } else {
            kmers_encoded[i] = kmer_to_uint32(kmers_list[i], k);
        }
    }

    // Compute Hamming distances
    for (size_t i = 0; i < n; i++) {
        uint32_t xor = query_kmer ^ kmers_encoded[i];
        uint32_t mask = (xor | (xor >> 1)) & 0x55555555u;

#if defined(__GNUC__) || defined(__clang__)
        distances[i] = __builtin_popcount(mask);
#else
        int count = 0;
        while (mask) {
            mask &= (mask - 1);
            count++;
        }
        distances[i] = count;
#endif
    }

    free(kmers_encoded);
}

void hamming_distance_array_64bit(
    const char *query,
    const char **kmers_list,
    size_t n,
    int *distances) 
{
    int k = strlen(query);

    // Encode query k-mer once
    uint64_t query_kmer = kmer_to_uint64(query, k);

    // Allocate temporary array for all kmers encoded
    uint64_t* kmers_encoded = (uint64_t*)malloc(n * sizeof(uint64_t));
    if (!kmers_encoded) return;  // allocation failure

    // Encode all kmers first
    for (size_t i = 0; i < n; i++) {
        if (!kmers_list[i]) {
            kmers_encoded[i] = 0;
        } else {
            kmers_encoded[i] = kmer_to_uint64(kmers_list[i], k);
        }
    }

    // Compute Hamming distances
    for (size_t i = 0; i < n; i++) {
        uint64_t xor = query_kmer ^ kmers_encoded[i];
        uint64_t mask = (xor | (xor >> 1)) & 0x5555555555555555u;

#if defined(__GNUC__) || defined(__clang__)
        distances[i] = __builtin_popcountll(mask);
#else
        int count = 0;
        while (mask) {
            mask &= (mask - 1);
            count++;
        }
        distances[i] = count;
#endif
    }

    free(kmers_encoded);
}

uint32_t hamming_distance_32bit(const char* kmer1_str, const char* kmer2_str) {
    int k = strlen(kmer1_str);  // assume kmer1_str and kmer2_str are same length
    uint32_t kmer1 = kmer_to_uint32(kmer1_str, k);
    uint32_t kmer2 = kmer_to_uint32(kmer2_str, k);

    uint32_t xor = kmer1 ^ kmer2;
    uint32_t mask = (xor | (xor >> 1)) & 0x55555555u;

#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(mask);
#else
    int count = 0;
    while (mask) {
        mask &= (mask - 1);
        count++;
    }
    return count;
#endif
}

uint64_t hamming_distance_64bit(const char* kmer1_str, const char* kmer2_str) {
    int k = strlen(kmer1_str);  // assume kmer1_str and kmer2_str are same length
    uint64_t kmer1 = kmer_to_uint64(kmer1_str, k);
    uint64_t kmer2 = kmer_to_uint64(kmer2_str, k);

    uint64_t xor = kmer1 ^ kmer2;
    uint64_t mask = (xor | (xor >> 1)) & 0x5555555555555555u;

#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(mask);
#else
    int count = 0;
    while (mask) {
        mask &= (mask - 1);
        count++;
    }
    return count;
#endif
}
