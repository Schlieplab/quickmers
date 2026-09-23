/*
Filename: quickmers/include/hamming.h
Author: Kian Jalilian
Copyright: 2025, Alexander Schliep
Version: 0.1.2
Description: Header file for hamming distance implementation in C
License: LGPL-3.0-or-later
*/
#ifndef HAMMING_H
#define HAMMING_H

#include <stdint.h>
#include <stddef.h>

void hamming_distance_array_32bit(const char *query, const char **kmers_list, size_t n, int *distances);

void hamming_distance_array_64bit(const char *query, const char **kmers_list, size_t n, int *distances);

uint64_t hamming_distance_64bit(const char* kmer1_str, const char* kmer2_str);

#endif