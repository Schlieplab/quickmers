#include "levenshtein.h"
#include <Python.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>

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

// 

// int64_t myers(uint8_t *t, int64_t n, uint8_t *p, int64_t m) {
//     uint64_t pv, mv;        // positive and negative vertical delta values
//     uint64_t ph, mh;        // positive and negative horizontal delta values
//     uint64_t xv, xh;        // current vertical and horizontal states
//     uint64_t eq, hb, peq[256];
//     int64_t score = m;
//     int64_t j;

//     // populate lookup table
//     memset(peq, 0, sizeof(peq));
//     for (int64_t i = 0; i < m; i++) {
//         // peq[p[i]] |= (1ULL << i);
//         peq[p[i]] |= (uint64_t) 1 << i;
//     }

//     // set initial values
//     pv = (1ULL << m) - 1;   // all ones up to pattern length
//     mv = 0ULL;
//     hb = 1ULL << (m - 1);   // high bit mask

//     // printf("pv after initialization: ");
//     // print_uint64_binary(pv);
//     // printf("mv after initialization: ");
//     // print_uint64_binary(mv);
//     // printf("hb after initialization: ");
//     // print_uint64_binary(hb);
//     // fflush(stdout); // ensures output appears immediately

//     for (j = 0; j < n; j++) {
//         eq = peq[t[j]];

//         // compute current states
//         xv = eq | mv;
//         xh = (((eq & pv) + pv) ^ pv) | eq;

//         // horizontal delta updates
//         ph = mv | ~(xh | pv);
//         mh = pv & xh;

//         // update score
//         if (ph & hb) {
//             score++;
//         } else if (mh & hb) {
//             score--;
//         }

//         // vertical delta updates
//         ph = (ph << 1) | 1ULL;
//         mh <<= 1;

//         pv = mh | ~(xv | ph);
//         mv = ph & xv;
//     }

//     return score;
// }