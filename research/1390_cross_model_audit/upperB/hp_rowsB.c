/* hp_rowsB.c -- 1390 (HP): float filter over all formal rows of Paper B, Proposition 9 (finitely many cases),
 * written from Paper B Sections 2-4 only.  Every move and every predecessor is made literally on a padded
 * word; W and Z are computed from the whole padded words (window table U times the gram factors V of the
 * left and right sides), and the ratio is formed after cancelling equal gram factors exactly by multiset
 * difference.  The row value is
 *     phi = 1/2 sum_{succ} theta(T->T') W(T')/W(T) + 1/2 sum_{pred} theta(T->T') Z(T')/Z(T).
 * Double arithmetic; every table value is correctly rounded from the exact rational (done in Python).
 * Rows with float value above (Rref - margin) are written out for exact re-evaluation in Python.
 * Internal codes: letter + = 1, - = 0; every word code is MSB = leftmost letter (grams: MSB = first letter
 * in reading order, i.e. from the far end toward the head).
 * Usage: hp_rowsB tables.bin Rref margin out_flagged.txt [threads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

static int r, k, q, d;
static double *Ut, *Ft, *Et, *Gt;
static double Rref, margin;

#define PADW 13
#define MAXL 64

typedef struct { uint8_t w[MAXL]; int len, h; } tape_t;

static uint32_t code_of(const uint8_t *a, int n) { uint32_t c = 0; for (int i = 0; i < n; i++) c = (c << 1) | a[i]; return c; }

static void grams(const tape_t *t, int comp, int sign, int16_t *cnt, int *touched, int *nt) {
    /* left side: positions 0..h-1, read left to right; right side: positions h+1..len-1, read right to left */
    for (int j = 0; j + k <= t->h; j++) {
        uint32_t c = 0; for (int i = 0; i < k; i++) c = (c << 1) | (uint32_t)(t->w[j + i] ^ comp);
        if (cnt[c] == 0) touched[(*nt)++] = (int)c;
        cnt[c] += sign;
    }
    for (int j = t->h + 1; j + k <= t->len; j++) {
        uint32_t c = 0; for (int i = 0; i < k; i++) c = (c << 1) | (uint32_t)(t->w[j + k - 1 - i] ^ comp);
        if (cnt[c] == 0) touched[(*nt)++] = (int)c;
        cnt[c] += sign;
    }
}
static uint32_t win_of(const tape_t *t, int comp, int rad) {
    uint32_t c = 0; for (int j = t->h - rad; j <= t->h + rad; j++) c = (c << 1) | (uint32_t)(t->w[j] ^ comp);
    return c;
}
/* ratio W(b)/W(a) (comp=0) or Z(b)/Z(a) (comp=1) */
static double ratio(const tape_t *a, const tape_t *b, int comp) {
    static __thread int16_t cnt[1 << 12];
    int touched[256], nt = 0;
    grams(b, comp, +1, cnt, touched, &nt);
    grams(a, comp, -1, cnt, touched, &nt);
    double x = Ut[win_of(b, comp, r)] / Ut[win_of(a, comp, r)];
    for (int i = 0; i < nt; i++) {
        int c = touched[i]; int e = cnt[c];
        /* the same code can be touched twice (count went to 0 and back); handle once */
        if (e > 0) for (int z = 0; z < e; z++) x *= Ft[c];
        else if (e < 0) for (int z = 0; z < -e; z++) x /= Ft[c];
        cnt[c] = 0;
    }
    return x;
}
static double cutE(const tape_t *t, int a0) { /* E of the 2q letters starting at position a0 */
    uint32_t c = 0; for (int i = 0; i < 2 * q; i++) c = (c << 1) | t->w[a0 + i];
    return Et[c];
}
static double Gof(const tape_t *shortt, int t0, int t1, int t2) {
    uint32_t wv = win_of(shortt, 0, d);
    return Gt[((uint64_t)wv << 3) | (uint32_t)(t0 << 2 | t1 << 1 | t2)];
}
static void replace1to3(const tape_t *a, tape_t *b, int x, int y, int z) { /* head letter -> x y z, head on y */
    int h = a->h; memcpy(b->w, a->w, h); b->w[h] = x; b->w[h + 1] = y; b->w[h + 2] = z;
    memcpy(b->w + h + 3, a->w + h + 1, a->len - h - 1); b->len = a->len + 2; b->h = h + 1;
}
static void replace3to1(const tape_t *a, tape_t *b, int x) { /* letters h-1,h,h+1 -> x, head on it */
    int h = a->h; memcpy(b->w, a->w, h - 1); b->w[h - 1] = x;
    memcpy(b->w + h, a->w + h + 2, a->len - h - 2); b->len = a->len - 2; b->h = h - 1;
}

static double row_value(const tape_t *T, int *nterms) {
    double s = 0; int n = 0; tape_t N; int h = T->h;
    /* successors */
    if (h >= 1) { N = *T; N.h = h - 1; s += (1.0 / cutE(T, h - q)) * ratio(T, &N, 0); n++; }          /* L: 1/E(cut left of head) */
    if (h + 1 < T->len) { N = *T; N.h = h + 1; s += cutE(T, h - q + 1) * ratio(T, &N, 0); n++; }     /* R: E(cut right of head) */
    if (T->w[h] == 1) { replace1to3(T, &N, 1, 0, 1); s += Gof(T, 1, 0, 1) * ratio(T, &N, 0); n++; } /* S1 */
    else {
        replace1to3(T, &N, 1, 0, 0); s += Gof(T, 1, 0, 0) * ratio(T, &N, 0); n++;                    /* S2 */
        replace1to3(T, &N, 0, 0, 1); s += Gof(T, 0, 0, 1) * ratio(T, &N, 0); n++;                    /* S3 */
    }
    if (h >= 1 && h + 1 < T->len && T->w[h] == 1 && !(T->w[h - 1] == 1 && T->w[h + 1] == 1)) {     /* merge */
        int x = T->w[h - 1], y = T->w[h + 1];
        int m = ((x ? 1 : -1) + 1 + (y ? 1 : -1)) > 0 ? 1 : 0;
        replace3to1(T, &N, m);
        s += (1.0 / Gof(&N, x, 1, y)) * ratio(T, &N, 0); n++;
    }
    /* predecessors */
    double p = 0;
    if (h + 1 < T->len) { N = *T; N.h = h + 1; p += cutE(T, h - q + 1) * ratio(T, &N, 1); n++; }    /* N moved L into T; theta(T->N) = right move */
    if (h >= 1) { N = *T; N.h = h - 1; p += (1.0 / cutE(T, h - q)) * ratio(T, &N, 1); n++; }         /* N moved R into T; theta = left move */
    if (h >= 1 && h + 1 < T->len && T->w[h] == 0) {                                                 /* split predecessors */
        int x = T->w[h - 1], y = T->w[h + 1], sl = -1;
        if (x == 1 && y == 1) sl = 1; else if (x == 1 && y == 0) sl = 0; else if (x == 0 && y == 1) sl = 0;
        if (sl >= 0) { replace3to1(T, &N, sl); p += (1.0 / Gof(&N, x, 0, y)) * ratio(T, &N, 1); n++; }
    }
    if (T->w[h] == 0) { replace1to3(T, &N, 0, 1, 0); p += Gof(T, 0, 1, 0) * ratio(T, &N, 1); n++; } /* merge predecessors */
    else {
        replace1to3(T, &N, 0, 1, 1); p += Gof(T, 0, 1, 1) * ratio(T, &N, 1); n++;
        replace1to3(T, &N, 1, 1, 0); p += Gof(T, 1, 1, 0) * ratio(T, &N, 1); n++;
    }
    *nterms = n;
    return 0.5 * s + 0.5 * p;
}

typedef struct { int tid, nth; double mx; uint32_t arg; uint64_t nflag; FILE *out; pthread_mutex_t *mu; uint64_t hist[8]; } job_t;

static void *worker(void *vp) {
    job_t *j = (job_t *)vp;
    int L = 2 * r + 3, half = r + 1;
    uint32_t nrows = 1u << L;
    for (uint32_t row = (uint32_t)j->tid; row < nrows; row += (uint32_t)j->nth) {
        tape_t T; memset(&T, 0, sizeof T);
        T.len = PADW + L + PADW; T.h = PADW + half;
        for (int i = 0; i < L; i++) T.w[PADW + i] = (uint8_t)((row >> (L - 1 - i)) & 1);
        int nt; double v = row_value(&T, &nt);
        if (v > j->mx) { j->mx = v; j->arg = row; }
        if (v > Rref - margin) {
            pthread_mutex_lock(j->mu); fprintf(j->out, "%u %.17g\n", row, v); pthread_mutex_unlock(j->mu); j->nflag++;
        }
        double gap = Rref - v; int b = gap < 0 ? 0 : gap < 1e-9 ? 1 : gap < 1e-6 ? 2 : gap < 1e-3 ? 3 : 4;
        j->hist[b]++;
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage\n"); return 2; }
    FILE *f = fopen(argv[1], "rb"); int32_t hd[4];
    if (fread(hd, 4, 4, f) != 4) return 2;
    r = hd[0]; k = hd[1]; q = hd[2]; d = hd[3];
    if (!(k <= r + 1 && q <= r + 1 && d <= r)) { printf("locality conditions fail\n"); return 1; }
    size_t nU = 1ull << (2 * r + 1), nF = 1ull << k, nE = 1ull << (2 * q), nG = (1ull << (2 * d + 1)) * 8;
    Ut = malloc(nU * 8); Ft = malloc(nF * 8); Et = malloc(nE * 8); Gt = malloc(nG * 8);
    if (fread(Ut, 8, nU, f) != nU || fread(Ft, 8, nF, f) != nF || fread(Et, 8, nE, f) != nE || fread(Gt, 8, nG, f) != nG) { printf("short tables\n"); return 2; }
    fclose(f);
    Rref = atof(argv[2]); margin = atof(argv[3]);
    int nth = argc > 5 ? atoi(argv[5]) : 4;
    FILE *out = fopen(argv[4], "w"); pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
    pthread_t th[16]; job_t jobs[16];
    for (int i = 0; i < nth; i++) { memset(&jobs[i], 0, sizeof(job_t)); jobs[i].tid = i; jobs[i].nth = nth; jobs[i].out = out; jobs[i].mu = &mu; pthread_create(&th[i], 0, worker, &jobs[i]); }
    double mx = 0; uint32_t arg = 0; uint64_t nflag = 0, hist[8] = {0};
    for (int i = 0; i < nth; i++) {
        pthread_join(th[i], 0); if (jobs[i].mx > mx) { mx = jobs[i].mx; arg = jobs[i].arg; } nflag += jobs[i].nflag;
        for (int b = 0; b < 8; b++) hist[b] += jobs[i].hist[b];
    }
    fclose(out);
    printf("r=%d k=%d q=%d d=%d rows=%u\n", r, k, q, d, 1u << (2 * r + 3));
    printf("float max row value %.15f at row %u (word ", mx, arg);
    for (int i = 0; i < 2 * r + 3; i++) { if (i == r + 1) putchar('['); putchar(((arg >> (2 * r + 2 - i)) & 1) ? '+' : '-'); if (i == r + 1) putchar(']'); }
    printf(")\nrows above Rref-margin: %llu\n", (unsigned long long)nflag);
    printf("gap histogram (Rref - v): <0:%llu [0,1e-9):%llu [1e-9,1e-6):%llu [1e-6,1e-3):%llu >=1e-3:%llu\n",
           (unsigned long long)hist[0], (unsigned long long)hist[1], (unsigned long long)hist[2], (unsigned long long)hist[3], (unsigned long long)hist[4]);
    return 0;
}
