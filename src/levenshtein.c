#include "levenshtein.h"
#include <Python.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <immintrin.h>  // AVX2

// print 64-bit integer as binary
void print_uint64_binary(uint64_t x) {
    for (int i = 63; i >= 0; i--) {
        putchar((x & ((uint64_t)1 << i)) ? '1' : '0');
    }
    putchar('\n');
}

/*
Levenshtein edit distance calculated via Myers bit-parallel algorithm

G. Myers (1999) A Fast Bit-Vector Algorithm for Approximate String Matching
Based on Dynamic Programming, Journal of the ACM: 495-415.
*/
int64_t myers(uint8_t *t, int64_t n, uint8_t *p, int64_t m) {
    uint64_t pv, mv;        // positive and negative vertical delta values
    uint64_t ph, mh;        // positive and negative horizontal delta values
    uint64_t xv, xh;        // current vertical and horizontal states
    uint64_t eq, hb, peq[256];
    int64_t score =  m;
    uint8_t i, j;

    // populate lookup table
    memset(peq, 0, sizeof(peq));
    for (i=0; i<m; i++) {
        peq[p[i]] |= (uint64_t) 1 << i;
    }

    // set initial values
    pv = (uint64_t) - 1;
    mv = (uint64_t) 0;
    hb = (uint64_t) 1 << (m - 1);
    
    for (j=0; j<n; j++) {
        // bit vectors for reference symbol j
        eq = peq[t[j]];

        // compute current delta vector
        xv = eq | mv;
        xh = (((eq & pv) + pv) ^ pv) | eq;

        // update horizontal delta values
        ph = mv | ~(xh | pv);
        mh = pv & xh;

        // update score
        if (ph & hb) {
            score++;
        }

        if (mh & hb) {
            score--;
        }

        // update vertical delta values
        ph = (ph << 1) | 1;
        pv = (mh << 1) | ~ (xv | ph);
        mv = ph & xv;
    }
    return score;
}

// AVX2 batch Myers: 4 kmers at once
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

        // Prepare peq tables for 4 kmers
        uint64_t peq[4][256] = {{0}};
        for (int lane = 0; lane < 4; lane++) {
            if (batch + lane >= n_kmers) break;
            const uint8_t *p = kmers[batch + lane];
            for (int i = 0; i < kmer_len; i++)
                peq[lane][p[i]] |= 1ULL << i;
        }

        // Loop over query positions
        for (int64_t j = 0; j < qlen; j++) {
            uint64_t eq_arr[4] = {0};
            for (int lane = 0; lane < 4; lane++) {
                if (batch + lane >= n_kmers) break;
                eq_arr[lane] = peq[lane][query[j]];
            }

            eq = _mm256_set_epi64x(
                eq_arr[3], eq_arr[2], eq_arr[1], eq_arr[0]
            );

            xv = _mm256_or_si256(eq, mv);
            xh = _mm256_or_si256(
                    _mm256_xor_si256(
                        _mm256_add_epi64(_mm256_and_si256(eq, pv), pv), pv
                    ),
                    eq
                );

            ph = _mm256_or_si256(mv, _mm256_andnot_si256(_mm256_or_si256(xh, pv), _mm256_set1_epi64x(-1LL)));
            mh = _mm256_and_si256(pv, xh);

            // Update scores per lane
            uint64_t ph_arr[4], mh_arr[4];
            _mm256_storeu_si256((__m256i*)ph_arr, ph);
            _mm256_storeu_si256((__m256i*)mh_arr, mh);
            for (int lane = 0; lane < 4; lane++) {
                if (batch + lane >= n_kmers) break;
                if (ph_arr[lane] & hb) score[lane]++;
                if (mh_arr[lane] & hb) score[lane]--;
            }

            // Update pv, mv
            ph = _mm256_slli_epi64(ph, 1);
            ph = _mm256_or_si256(ph, _mm256_set1_epi64x(1));
            pv = _mm256_or_si256(
                    _mm256_slli_epi64(mh,1),
                    _mm256_andnot_si256(_mm256_or_si256(xv, ph), _mm256_set1_epi64x(-1LL))
                 );
            mv = _mm256_and_si256(ph, xv);
        }

        // Store results
        for (int lane = 0; lane < 4; lane++) {
            if (batch + lane >= n_kmers) break;
            out[batch + lane] = score[lane];
        }
    }
}