/*
Filename: src/levenshtein.c
Author: Kian Jalilian
Copyright: 2025, Alexander Schliep
Version: 0.2.0
Description: Function implementations for levenshtein distance calculations
License: LGPL-3.0-or-later
*/
#include "levenshtein.h"
#include <Python.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>

#if defined(__x86_64__) && defined(__AVX2__)
#include <immintrin.h>  // AVX2
#include <cpuid.h>
#endif

// ---------------------- Scalar Myers -------------------------

int64_t myers(uint8_t *t, int64_t n, uint8_t *p, int64_t m) {
    uint64_t pv = (uint64_t)-1;
    uint64_t mv = 0;
    uint64_t eq, hb, ph, mh, xv, xh, peq[256];
    int64_t score = m;
    uint8_t i, j;

    memset(peq, 0, sizeof(peq));
    for (i = 0; i < m; i++) peq[p[i]] |= (uint64_t)1 << i;

    hb = (uint64_t)1 << (m - 1);

    for (j = 0; j < n; j++) {
        eq = peq[t[j]];
        xv = eq | mv;
        xh = (((eq & pv) + pv) ^ pv) | eq;
        ph = mv | ~(xh | pv);
        mh = pv & xh;
        if (ph & hb) score++;
        if (mh & hb) score--;
        ph = (ph << 1) | 1;
        pv = (mh << 1) | ~(xv | ph);
        mv = ph & xv;
    }
    return score;
}

// ---------------------- AVX2 Myers ---------------------------

#if defined(__x86_64__) && defined(__AVX2__)
void myers_batch_avx2(
    const uint8_t *query, int64_t qlen,
    const uint8_t **kmers, int64_t n_kmers, int64_t kmer_len,
    int64_t *out
) {
    for (int64_t batch = 0; batch < n_kmers; batch += 4) {
        __m256i pv = _mm256_set1_epi64x(-1LL);
        __m256i mv = _mm256_setzero_si256();
        __m256i eq, xv, xh, ph, mh;
        int64_t score[4] = {kmer_len, kmer_len, kmer_len, kmer_len};
        uint64_t hb = 1ULL << (kmer_len - 1);

        uint64_t peq[4][256] = {{0}};
        for (int lane = 0; lane < 4; lane++) {
            if (batch + lane >= n_kmers) break;
            const uint8_t *p = kmers[batch + lane];
            for (int i = 0; i < kmer_len; i++) peq[lane][p[i]] |= 1ULL << i;
        }

        for (int64_t j = 0; j < qlen; j++) {
            uint64_t eq_arr[4] = {0};
            for (int lane = 0; lane < 4; lane++)
                if (batch + lane < n_kmers) eq_arr[lane] = peq[lane][query[j]];

            eq = _mm256_set_epi64x(eq_arr[3], eq_arr[2], eq_arr[1], eq_arr[0]);
            xv = _mm256_or_si256(eq, mv);
            xh = _mm256_or_si256(_mm256_xor_si256(_mm256_add_epi64(_mm256_and_si256(eq, pv), pv), pv), eq);

            ph = _mm256_or_si256(mv, _mm256_andnot_si256(_mm256_or_si256(xh, pv), _mm256_set1_epi64x(-1LL)));
            mh = _mm256_and_si256(pv, xh);

            uint64_t ph_arr[4], mh_arr[4];
            _mm256_storeu_si256((__m256i*)ph_arr, ph);
            _mm256_storeu_si256((__m256i*)mh_arr, mh);
            for (int lane = 0; lane < 4; lane++) {
                if (batch + lane >= n_kmers) break;
                if (ph_arr[lane] & hb) score[lane]++;
                if (mh_arr[lane] & hb) score[lane]--;
            }

            ph = _mm256_or_si256(_mm256_slli_epi64(ph, 1), _mm256_set1_epi64x(1));
            pv = _mm256_or_si256(_mm256_slli_epi64(mh, 1),
                                  _mm256_andnot_si256(_mm256_or_si256(xv, ph), _mm256_set1_epi64x(-1LL)));
            mv = _mm256_and_si256(ph, xv);
        }

        for (int lane = 0; lane < 4; lane++)
            if (batch + lane < n_kmers) out[batch + lane] = score[lane];
    }
}
#endif

// -------------------- AVX2 Detection ------------------------

static int has_avx2() {
#if defined(__x86_64__) && defined(__AVX2__)
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        return (ebx & bit_AVX2) != 0;
    }
#endif
    return 0;
}

// -------------------- Dispatcher -----------------------------

void myers_dispatch(
    const uint8_t *query, int64_t qlen,
    const uint8_t **kmers, int64_t n_kmers, int64_t kmer_len,
    int64_t *out
) {
    if (has_avx2()) {
#ifdef __x86_64__
        myers_batch_avx2(query, qlen, kmers, n_kmers, kmer_len, out);
        return;
#endif
    }
    // fallback scalar implementation
    for (int64_t i = 0; i < n_kmers; i++) {
        out[i] = myers((uint8_t*)kmers[i], kmer_len, (uint8_t*)query, qlen);
    }
}
