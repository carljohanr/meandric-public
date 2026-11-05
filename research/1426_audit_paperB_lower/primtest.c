/* 1426 (JB): unit test of the C primitives of flowcheck.c (own logarithm, beta mantissas, beta mod p).
 * Reads lines from stdin:  "L A_hi A_lo B_hi B_lo"  (two u128 as 64-bit halves) -> prints log_ratio_lo as decimal Q64
 *                          "B lev code len"           -> prints beta lo m e, hi m e, mod0, mod1
 * Checked against exact values by primtest.py. */
#define NOMAIN
#include "flowcheck.c"
int main(int argc, char **argv) {
    load_law(argc > 1 ? argv[1] : "work/law_c.txt");
    char t[4];
    while (scanf("%3s", t) == 1) {
        if (t[0] == 'L') {
            unsigned long long ah, al, bh, bl; if (scanf("%llu %llu %llu %llu", &ah, &al, &bh, &bl) != 4) break;
            u128 A = ((u128)ah << 64) | al, B = ((u128)bh << 64) | bl;
            u128 r = log_ratio_lo(A, B);
            print_u128("L", r);
        } else {
            int lev, len; unsigned long long c; if (scanf("%d %llu %d", &lev, &c, &len) != 3) break;
            mant lo, hi; beta_mant(lev, c, len, &lo, &hi); uint64_t m[2]; beta_mod(lev, c, len, m);
            printf("B %llu %d %llu %d %llu %llu\n", (unsigned long long)lo.m, lo.e, (unsigned long long)hi.m, hi.e,
                   (unsigned long long)m[0], (unsigned long long)m[1]);
        }
    }
    return 0;
}
