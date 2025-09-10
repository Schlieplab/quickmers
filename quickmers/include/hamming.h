#ifndef HAMMING_H
#define HAMMING_H

#include <stdint.h>
#include <stddef.h>

void hamming_distance_encoded_array_32bit(uint32_t kmer, const uint32_t *kmers_list, size_t n, int *distances);
void hamming_distance_array_32bit(const char *query, const char **kmers_list, size_t n, int *distances);

void hamming_distance_encoded_array_64bit(uint64_t kmer, const uint64_t *kmers_list, size_t n, int *distances);
void hamming_distance_array_64bit(const char *query, const char **kmers_list, size_t n, int *distances);

uint32_t hamming_distance_encoded_32bit(uint32_t kmer1, uint32_t kmer2);
uint64_t hamming_distance_encoded_64bit(uint64_t kmer1, uint64_t kmer2);

#endif