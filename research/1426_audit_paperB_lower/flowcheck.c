/* 1426 (JB): independent checker for the final flow of Paper B's lower bound (SPEC 1421 Sections 5-7).
 *
 * Builds J = J_P lifted + sum_d c_d Gamma_d at an evaluation radius R from
 *   work/law_c.txt    (law tables: quartets q, quintet factors t = p/q; mod p in Montgomery form; 64-bit mantissas)
 *   work/parent_c.bin (parent densities x_P at radius 6: floor(x*2^88) and x mod p)
 *   work/edges_c.bin  (one record per walk edge: source tape, move, stage, term c*beta(A)/beta(T))
 * Densities are kept as exact fixed-point sums: each term is floor(term*2^88), so for a pair hit by n terms
 *   lo <= x*2^88 <= lo + n.
 * Checks, at radius R, on EVERY legal pair / balance row:
 *   (P)  x_lo > 0 on all legal pairs;
 *   (S)  every balance row (S_R) of SPEC 7.1 multiplied by beta(Y), modulo two random 62-bit primes (exact in Z/p);
 *   (D)  sum beta*x*Da = sum beta*x*Db = 0 mod p; mass sum beta*x mod p (compared with the closed form by bound.py);
 *   entropy: H_lo = sum beta_lo x_lo log_lo(nubar_lo / x_hi) and M_hi = sum beta_hi nubar_hi, in Q120 fixed point,
 *   with directed rounding everywhere and an own logarithm (atanh series, rigorous tail).
 * Usage: flowcheck R SMAX THREADS ROWS(0/1) [parent file] [edges file] [law file]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

typedef unsigned __int128 u128;
typedef __int128 i128;
#define K 12
#define SCALE 88

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void die(const char *m) { fprintf(stderr, "FATAL: %s\n", m); exit(2); }

/* ---------------- Montgomery arithmetic, two primes ---------------- */
static uint64_t PR[2], PINV[2]; /* PINV = -p^{-1} mod 2^64 */
static inline uint64_t mmul(uint64_t a, uint64_t b, int k) {
    u128 T = (u128)a * b;
    uint64_t m = (uint64_t)T * PINV[k];
    u128 t = (T + (u128)m * PR[k]) >> 64;
    uint64_t r = (uint64_t)t;
    if (r >= PR[k]) r -= PR[k];
    return r;
}
static inline uint64_t madd(uint64_t a, uint64_t b, int k) { uint64_t r = a + b; if (r >= PR[k]) r -= PR[k]; return r; }
static inline uint64_t msub(uint64_t a, uint64_t b, int k) { return a >= b ? a - b : a + PR[k] - b; }

/* ---------------- 64-bit mantissa numbers with directed rounding ---------------- */
typedef struct { uint64_t m; int e; } mant; /* value m*2^e, m in [2^63,2^64) */
static inline mant mul_dn(mant a, mant b) {
    u128 P = (u128)a.m * b.m; mant r;
    if (P >> 127) { r.m = (uint64_t)(P >> 64); r.e = a.e + b.e + 64; }
    else { r.m = (uint64_t)(P >> 63); r.e = a.e + b.e + 63; }
    return r;
}
static inline mant mul_up(mant a, mant b) {
    u128 P = (u128)a.m * b.m; mant r; int s;
    if (P >> 127) s = 64; else s = 63;
    u128 q = P >> s; int rem = (P & ((((u128)1) << s) - 1)) != 0;
    q += rem;
    r.e = a.e + b.e + s;
    if (q >> 64) { r.m = (uint64_t)1 << 63; r.e += 1; } else r.m = (uint64_t)q;
    return r;
}
static inline int bitlen128(u128 x) {
    uint64_t hi = (uint64_t)(x >> 64), lo = (uint64_t)x;
    if (hi) return 128 - __builtin_clzll(hi);
    if (lo) return 64 - __builtin_clzll(lo);
    return 0;
}
/* x > 0, value x*2^ex */
static inline mant from_u128_dn(u128 x, int ex) {
    int bl = bitlen128(x); mant r;
    if (bl > 64) { r.m = (uint64_t)(x >> (bl - 64)); r.e = ex + bl - 64; }
    else { r.m = (uint64_t)x << (64 - bl); r.e = ex - (64 - bl); }
    return r;
}
static inline mant from_u128_up(u128 x, int ex) {
    int bl = bitlen128(x); mant r;
    if (bl > 64) {
        int s = bl - 64; u128 q = x >> s; int rem = (x & ((((u128)1) << s) - 1)) != 0; q += rem;
        r.e = ex + s;
        if (q >> 64) { r.m = (uint64_t)1 << 63; r.e += 1; } else r.m = (uint64_t)q;
    } else { r.m = (uint64_t)x << (64 - bl); r.e = ex - (64 - bl); }
    return r;
}
/* floor / ceil of value * 2^120 as u128; requires value < 2^7 */
static inline u128 to_q120_dn(mant a) {
    int s = a.e + 120;
    if (s >= 0) { if (s > 63) die("q120 overflow"); return ((u128)a.m) << s; }
    if (-s >= 128) return 0;
    return ((u128)a.m) >> (-s);
}
static inline u128 to_q120_up(mant a) {
    int s = a.e + 120;
    if (s >= 0) { if (s > 63) die("q120 overflow"); return ((u128)a.m) << s; }
    if (-s >= 128) return 1;
    u128 q = ((u128)a.m) >> (-s);
    if (((u128)a.m) & ((((u128)1) << (-s)) - 1)) q += 1;
    return q;
}

/* ---------------- own logarithm: lower bound of ln(A/B) in Q64 ---------------- */
static const uint64_t LN2_LO = 0xb17217f7d1cf79abULL; /* floor(ln2 * 2^64), checked two ways (logs/constants) */
static const uint64_t LN2_HI = 0xb17217f7d1cf79acULL;
/* atanh(u) for u = U/2^64 in [0, 0.2]: lower (dir=0) or upper (dir=1) bound, Q64 */
static inline u128 atanh_q64(u128 U, int dir) {
    if (U == 0) return 0;
    u128 u2 = U * U; /* U < 2^63 so U*U < 2^126 fits */
    u128 sq = u2 >> 64; if (dir && (u2 & (u128)0xFFFFFFFFFFFFFFFFULL)) sq += 1;
    u128 pw = U, sum = 0;
    for (int j = 0; j < 16; j++) {
        u128 d = 2 * j + 1;
        u128 t = pw / d; if (dir && (pw % d)) t += 1;
        sum += t;
        u128 nx = pw * sq; /* pw < 2^63, sq < 2^62 */
        u128 np = nx >> 64; if (dir && (nx & (u128)0xFFFFFFFFFFFFFFFFULL)) np += 1;
        pw = np;
    }
    if (dir) sum += 2 * pw + 1; /* tail: sum_{j>=16} u^{2j+1}/(2j+1) <= u^33/(1-u^2) <= 2 u^33 */
    return sum;
}
/* lower bound of ln(A/B) as Q64 (>= 0), A, B > 0 integers */
static inline u128 log_ratio_lo(u128 A, u128 B) {
    if (A <= B) return 0;
    int la = bitlen128(A), lb = bitlen128(B);
    u128 a, b; int ea, eb;
    if (la > 63) { a = A >> (la - 63); ea = la - 63; } else { a = A << (63 - la); ea = -(63 - la); }
    if (lb > 63) {
        int s = lb - 63; b = B >> s; if (B & ((((u128)1) << s) - 1)) b += 1; eb = s;
    } else { b = B << (63 - lb); eb = -(63 - lb); }
    /* A/B >= (a/b) 2^(ea-eb), a in [2^62,2^63), b in [2^62,2^63] */
    int e = ea - eb;
    double md = (double)a / (double)b;
    if (md < 0.7071) { a <<= 1; e -= 1; } else if (md > 1.4142) { b <<= 1; e += 1; }
    i128 res;
    if (a >= b) {
        u128 num = (a - b) << 64, den = a + b;
        u128 U = num / den;                        /* rounded down */
        if (U > (u128)3689348814741910323ULL) die("atanh argument above 0.2");
        res = (i128)(2 * atanh_q64(U, 0));
    } else {
        u128 num = (b - a) << 64, den = a + b;
        u128 U = num / den; if (num % den) U += 1;  /* rounded up */
        if (U > (u128)3689348814741910323ULL) die("atanh argument above 0.2");
        res = -(i128)(2 * atanh_q64(U, 1));
    }
    if (e >= 0) res += (i128)e * (i128)LN2_LO; else res += (i128)e * (i128)LN2_HI;
    if (res < 0) res = 0;
    return (u128)res;
}

/* ---------------- law tables ---------------- */
static uint64_t QMOD[2][96], TMOD[2][192];
static mant QLO[96], QHI[96], TLO[192], THI[192];
static int QOK[96], TOK[192];
static void load_law(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) die("law file");
    unsigned long long p0, p1;
    if (fscanf(f, "%llu %llu", &p0, &p1) != 2) die("law header");
    PR[0] = p0; PR[1] = p1;
    for (int k = 0; k < 2; k++) { uint64_t inv = 1; for (int i = 0; i < 7; i++) inv *= 2 - PR[k] * inv; PINV[k] = -inv; if (PR[k] * inv != 1) die("pinv"); }
    char t[4]; int idx; unsigned long long m0, m1, lo, hi; int elo, ehi; int nq = 0, nt = 0;
    while (fscanf(f, "%3s %d %llu %llu %llu %d %llu %d", t, &idx, &m0, &m1, &lo, &elo, &hi, &ehi) == 8) {
        if (t[0] == 'Q') { QMOD[0][idx] = m0; QMOD[1][idx] = m1; QLO[idx].m = lo; QLO[idx].e = elo; QHI[idx].m = hi; QHI[idx].e = ehi; QOK[idx] = 1; nq++; }
        else { TMOD[0][idx] = m0; TMOD[1][idx] = m1; TLO[idx].m = lo; TLO[idx].e = elo; THI[idx].m = hi; THI[idx].e = ehi; TOK[idx] = 1; nt++; }
        if ((lo >> 63) == 0 || (hi >> 63) == 0) die("mantissa not normalized");
    }
    fclose(f);
    if (nq != 80 || nt != 154) die("law table counts");
}
/* word: level, code (first letter = most significant bit, 1 = '+'), len letters (len >= 3) */
static inline void beta_mod(int lev, uint64_t c, int len, uint64_t out[2]) {
    int qi = lev * 8 + (int)((c >> (len - 3)) & 7);
    uint64_t r0 = QMOD[0][qi], r1 = QMOD[1][qi];
    int h = lev;
    for (int i = 4; i <= len; i++) {
        int ti = h * 16 + (int)((c >> (len - i)) & 15);
        r0 = mmul(r0, TMOD[0][ti], 0); r1 = mmul(r1, TMOD[1][ti], 1);
        h += ((c >> (len - (i - 3))) & 1) ? 1 : -1; /* height after letter i-4 */
    }
    out[0] = r0; out[1] = r1;
}
static inline void beta_mant(int lev, uint64_t c, int len, mant *lo, mant *hi) {
    int qi = lev * 8 + (int)((c >> (len - 3)) & 7);
    mant a = QLO[qi], b = QHI[qi];
    int h = lev;
    for (int i = 4; i <= len; i++) {
        int ti = h * 16 + (int)((c >> (len - i)) & 15);
        if (!TOK[ti]) die("quintet missing");
        a = mul_dn(a, TLO[ti]); b = mul_up(b, THI[ti]);
        h += ((c >> (len - (i - 3))) & 1) ? 1 : -1;
    }
    *lo = a; *hi = b;
}
static inline int in_strip(int lev, uint64_t c, int len) {
    int h = lev;
    if (h < 0 || h >= K) return 0;
    for (int i = len - 1; i >= 0; i--) { h += ((c >> i) & 1) ? 1 : -1; if (h < 0 || h >= K) return 0; }
    return 1;
}

/* ---------------- windows at a radius ---------------- */
typedef struct {
    int r, n; int32_t *widx; uint32_t nwin; uint32_t *wkey; uint32_t *wbase; uint8_t *wmask; uint64_t npairs;
    i128 *lo; uint32_t *cnt; uint64_t *m0, *m1; int has; uint8_t *touched;
} Rad;
static Rad RD[11];
static uint8_t legal_mask(int lev, uint64_t c, int n) {
    int r = (n - 1) / 2;
    int h = lev; for (int i = 0; i < r; i++) h += ((c >> (n - 1 - i)) & 1) ? 1 : -1; /* height before head */
    int head = (c >> (n - 1 - r)) & 1, lft = (c >> (n - r)) & 1, rgt = (c >> (n - 2 - r)) & 1;
    uint8_t m = 3;
    if (head) { m |= 4; if (!(lft && rgt)) m |= 32; }
    else { if (h + 1 <= K - 1) m |= 8; if (h - 2 >= 0) m |= 16; }
    return m;
}
static void build_rad(int r) {
    Rad *R = &RD[r]; R->r = r; R->n = 2 * r + 1; int n = R->n;
    uint64_t space = (uint64_t)K << n;
    R->widx = malloc(space * sizeof(int32_t)); if (!R->widx) die("alloc widx");
    uint32_t cntw = 0;
    for (uint64_t key = 0; key < space; key++) {
        int lev = (int)(key >> n); uint64_t c = key & ((1ULL << n) - 1);
        R->widx[key] = in_strip(lev, c, n) ? (int32_t)cntw++ : -1;
    }
    R->nwin = cntw;
    R->wkey = malloc(cntw * 4ULL); R->wbase = malloc(cntw * 4ULL); R->wmask = malloc(cntw);
    uint64_t base = 0;
    for (uint64_t key = 0; key < space; key++) {
        int32_t i = R->widx[key]; if (i < 0) continue;
        R->wkey[i] = (uint32_t)key; R->wbase[i] = (uint32_t)base;
        R->wmask[i] = legal_mask((int)(key >> n), key & ((1ULL << n) - 1), n);
        base += __builtin_popcount(R->wmask[i]);
    }
    R->npairs = base;
    R->lo = calloc(base, sizeof(i128)); R->cnt = calloc(base, 4); R->m0 = calloc(base, 8); R->m1 = calloc(base, 8);
    if (!R->lo || !R->cnt || !R->m0 || !R->m1) die("alloc comps");
    R->has = 1;
    printf("radius %d: windows %u legal pairs %llu\n", r, cntw, (unsigned long long)base); fflush(stdout);
}
static inline int64_t pair_index(Rad *R, int lev, uint64_t c, int mv) {
    int32_t i = R->widx[((uint64_t)lev << R->n) | c];
    if (i < 0) return -1;
    uint8_t m = R->wmask[i];
    if (!((m >> mv) & 1)) return -2;
    return (int64_t)R->wbase[i] + __builtin_popcount(m & ((1u << mv) - 1));
}

/* ---------------- globals for the passes ---------------- */
static int EVR, NTH, SMAX;
static const int STAGE_R[5] = {8, 9, 10, 10, 10};
static int USED[11]; /* radii holding components (r < EVR) */

typedef struct {
    int tid;
    /* pass 1 outputs */
    uint64_t nonpos, pairs, wins; double minx; uint32_t maxcnt;
    u128 H, Mhi, Mlo; uint64_t mass[2], da[2], db[2];
    /* pass 2 outputs */
    uint64_t rows, badrows[2], anomalies; uint64_t firstbad_key;
} Th;

/* pass 1: lift lower-radius components into radius EVR, then positivity, entropy, mass, drift */
static void *pass1(void *arg) {
    Th *T = arg; Rad *R = &RD[EVR]; int n = R->n;
    T->minx = 1e300; T->firstbad_key = ~0ULL;
    uint64_t w0 = (uint64_t)R->nwin * T->tid / NTH, w1 = (uint64_t)R->nwin * (T->tid + 1) / NTH;
    for (uint64_t wi = w0; wi < w1; wi++) {
        uint32_t key = R->wkey[wi]; int lev = key >> n; uint64_t c = key & ((1ULL << n) - 1);
        uint8_t mask = R->wmask[wi]; uint64_t base = R->wbase[wi];
        /* centers */
        int64_t cidx[11][6]; int rs[4], nr = 0;
        for (int r = 6; r < EVR; r++) if (USED[r]) {
            int d = EVR - r; uint64_t cc = (c >> d) & ((1ULL << (2 * r + 1)) - 1);
            int pc = __builtin_popcountll(c >> (n - d)); int cl = lev + 2 * pc - d;
            for (int mv = 0; mv < 6; mv++) if ((mask >> mv) & 1) {
                cidx[r][mv] = pair_index(&RD[r], cl, cc, mv);
                if (cidx[r][mv] < 0) die("center pair missing");
            }
            rs[nr++] = r;
        }
        i128 xlo[6]; u128 xhi[6]; int k = 0; int mvs[6];
        uint64_t xm0[6], xm1[6];
        for (int mv = 0; mv < 6; mv++) if ((mask >> mv) & 1) {
            uint64_t pi = base + k;
            for (int j = 0; j < nr; j++) {
                Rad *C = &RD[rs[j]]; int64_t q = cidx[rs[j]][mv];
                R->lo[pi] += C->lo[q]; R->cnt[pi] += C->cnt[q];
                R->m0[pi] = madd(R->m0[pi], C->m0[q], 0); R->m1[pi] = madd(R->m1[pi], C->m1[q], 1);
            }
            if (R->cnt[pi] > T->maxcnt) T->maxcnt = R->cnt[pi];
            xlo[k] = R->lo[pi]; mvs[k] = mv; xm0[k] = R->m0[pi]; xm1[k] = R->m1[pi];
            k++;
        }
        T->wins++;
        int bad = 0;
        for (int j = 0; j < k; j++) {
            T->pairs++;
            if (xlo[j] <= 0) { T->nonpos++; bad = 1; }
            double xd = (double)xlo[j] / 309485009821345068724781056.0; /* 2^88 */
            if (xd < T->minx) T->minx = xd;
        }
        mant blo, bhi; beta_mant(lev, c, n, &blo, &bhi);
        uint64_t bm[2]; beta_mod(lev, c, n, bm);
        /* mod p: mass and drift */
        uint64_t s0 = 0, s1 = 0;
        for (int j = 0; j < k; j++) {
            uint64_t t0 = mmul(bm[0], xm0[j], 0), t1 = mmul(bm[1], xm1[j], 1);
            s0 = madd(s0, t0, 0); s1 = madd(s1, t1, 1);
            int mv = mvs[j]; int dA = (mv == 0 || mv == 5) ? -1 : 1; int dB = (mv == 0 || mv >= 2) ? 1 : -1;
            if (mv == 5) dB = -1;
            /* Delta(L)=(-1,+1) Delta(R)=(+1,-1) splits (+1,+1) M (-1,-1) */
            T->da[0] = dA > 0 ? madd(T->da[0], t0, 0) : msub(T->da[0], t0, 0);
            T->da[1] = dA > 0 ? madd(T->da[1], t1, 1) : msub(T->da[1], t1, 1);
            T->db[0] = dB > 0 ? madd(T->db[0], t0, 0) : msub(T->db[0], t0, 0);
            T->db[1] = dB > 0 ? madd(T->db[1], t1, 1) : msub(T->db[1], t1, 1);
        }
        T->mass[0] = madd(T->mass[0], s0, 0); T->mass[1] = madd(T->mass[1], s1, 1);
        if (bad) continue; /* entropy only meaningful when positive; the run fails anyway */
        for (int j = 0; j < k; j++) if (((u128)xlo[j]) >> 120) die("density above 2^32");
        u128 nulo = 0, nuhi = 0;
        for (int j = 0; j < k; j++) { nulo += (u128)xlo[j]; nuhi += (u128)xlo[j] + R->cnt[base + j]; }
        for (int j = 0; j < k; j++) {
            u128 xh = (u128)xlo[j] + R->cnt[base + j];
            u128 Lq = log_ratio_lo(nulo, xh);
            if (Lq == 0) continue;
            mant xm = from_u128_dn((u128)xlo[j], -SCALE);
            mant lm = from_u128_dn(Lq, -64);
            mant t = mul_dn(mul_dn(blo, xm), lm);
            T->H += to_q120_dn(t);
        }
        T->Mhi += to_q120_up(mul_up(bhi, from_u128_up(nuhi, -SCALE)));
        T->Mlo += to_q120_dn(mul_dn(blo, from_u128_dn(nulo, -SCALE)));
    }
    return NULL;
}

/* pass 2: every balance row (S_EVR), mod both primes, multiplied by beta(Y) */
static void *pass2(void *arg) {
    Th *T = arg; Rad *R = &RD[EVR]; int r = EVR, n = R->n, m = n + 2;
    uint64_t space = (uint64_t)K << m;
    uint64_t k0 = space * T->tid / NTH, k1 = space * (T->tid + 1) / NTH;
    uint64_t maskn = (1ULL << n) - 1;
    for (uint64_t key = k0; key < k1; key++) {
        int lev = (int)(key >> m); uint64_t cY = key & ((1ULL << m) - 1);
        if (!in_strip(lev, cY, m)) continue;
        T->rows++;
        int y0 = (cY >> (m - 1)) & 1, y1 = (cY >> (m - 2)) & 1;
        int lev1 = lev + (y0 ? 1 : -1), lev2 = lev1 + (y1 ? 1 : -1);
        uint64_t bY[2]; beta_mod(lev, cY, m, bY);
        uint64_t acc0 = 0, acc1 = 0; /* x_L + x_R - sum x_u, times beta(Y) */
        int64_t q;
        q = pair_index(R, lev2, cY & maskn, 0); if (q < 0) { T->anomalies++; continue; }
        acc0 = madd(acc0, R->m0[q], 0); acc1 = madd(acc1, R->m1[q], 1);
        q = pair_index(R, lev, cY >> 2, 1); if (q < 0) { T->anomalies++; continue; }
        acc0 = madd(acc0, R->m0[q], 0); acc1 = madd(acc1, R->m1[q], 1);
        uint64_t cu = (cY >> 1) & maskn;
        int32_t ui = R->widx[((uint64_t)lev1 << n) | cu]; if (ui < 0) { T->anomalies++; continue; }
        { uint8_t um = R->wmask[ui]; uint64_t b = R->wbase[ui]; int kk = __builtin_popcount(um);
          for (int j = 0; j < kk; j++) { acc0 = msub(acc0, R->m0[b + j], 0); acc1 = msub(acc1, R->m1[b + j], 1); } }
        uint64_t res0 = mmul(acc0, bY[0], 0), res1 = mmul(acc1, bY[1], 1);
        /* split sources */
        int trip = (int)((cY >> r) & 7); int smv = -1, rep = 0;
        if (trip == 5) { smv = 2; rep = 1; } else if (trip == 4) { smv = 3; rep = 0; } else if (trip == 1) { smv = 4; rep = 0; }
        if (smv >= 0) {
            uint64_t cS = ((cY >> (r + 3)) << (r + 1)) | ((uint64_t)rep << r) | (cY & ((1ULL << r) - 1));
            q = pair_index(R, lev, cS, smv);
            if (q >= 0) {
                uint64_t bS[2]; beta_mod(lev, cS, n, bS);
                res0 = madd(res0, mmul(R->m0[q], bS[0], 0), 0); res1 = madd(res1, mmul(R->m1[q], bS[1], 1), 1);
            } else T->anomalies++; /* SPEC: term only if w_S is a window and S legal; never expected */
        }
        /* merge sources */
        int ym = (cY >> (r + 1)) & 1;
        for (int pq = 0; pq < 2; pq++) {
            int p, qq;
            if (!ym) { if (pq) break; p = 0; qq = 0; }
            else { p = pq ? 1 : 0; qq = pq ? 0 : 1; }
            uint64_t cV = ((cY >> (r + 2)) << (r + 4)) | ((uint64_t)p << (r + 3)) | (1ULL << (r + 2)) |
                          ((uint64_t)qq << (r + 1)) | (cY & ((1ULL << (r + 1)) - 1));
            if (!in_strip(lev, cV, m + 2)) continue;
            uint64_t cM = (cV >> 2) & maskn;
            q = pair_index(R, lev2, cM, 5);
            if (q < 0) { T->anomalies++; continue; }
            uint64_t bV[2]; beta_mod(lev, cV, m + 2, bV);
            res0 = madd(res0, mmul(R->m0[q], bV[0], 0), 0); res1 = madd(res1, mmul(R->m1[q], bV[1], 1), 1);
        }
        if (res0) { T->badrows[0]++; if (key < T->firstbad_key) T->firstbad_key = key; }
        if (res1) T->badrows[1]++;
    }
    return NULL;
}

static void print_u128(const char *name, u128 v) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (v == 0) buf[--i] = '0';
    while (v) { buf[--i] = '0' + (int)(v % 10); v /= 10; }
    printf("%s %s\n", name, buf + i);
}

#ifndef NOMAIN
int main(int argc, char **argv) {
    if (argc < 5) die("usage: flowcheck R SMAX THREADS ROWS [parent] [edges] [law]");
    EVR = atoi(argv[1]); SMAX = atoi(argv[2]); NTH = atoi(argv[3]); int ROWS = atoi(argv[4]);
    const char *fpar = argc > 5 ? argv[5] : "work/parent_c.bin";
    const char *fedg = argc > 6 ? argv[6] : "work/edges_c.bin";
    const char *flaw = argc > 7 ? argv[7] : "work/law_c.txt";
    double t0 = now();
    printf("flowcheck R=%d SMAX=%d threads=%d rows=%d parent=%s edges=%s law=%s\n", EVR, SMAX, NTH, ROWS, fpar, fedg, flaw);
    load_law(flaw);
    printf("primes %llu %llu\n", (unsigned long long)PR[0], (unsigned long long)PR[1]);
    /* radii with components */
    USED[6] = 1;
    for (int s = 0; s <= SMAX; s++) { if (STAGE_R[s] > EVR) die("stage radius above evaluation radius"); USED[STAGE_R[s]] = 1; }
    for (int r = 6; r <= EVR; r++) if (USED[r] || r == EVR) build_rad(r);
    printf("windows built %.1fs\n", now() - t0); fflush(stdout);
    /* parent */
    FILE *f = fopen(fpar, "rb"); if (!f) die("parent file");
    uint32_t np; if (fread(&np, 4, 1, f) != 1) die("parent header");
    double maxpar = 0;
    for (uint32_t i = 0; i < np; i++) {
        uint8_t hd[8]; i128 lo; uint64_t mm[2];
        if (fread(hd, 8, 1, f) != 1 || fread(&lo, 16, 1, f) != 1 || fread(mm, 16, 1, f) != 1) die("parent read");
        int lev = hd[0], mv = hd[1]; uint32_t cd; memcpy(&cd, hd + 4, 4);
        int64_t q = pair_index(&RD[6], lev, cd, mv); if (q < 0) die("parent pair not legal");
        if (RD[6].cnt[q]) die("parent pair twice");
        RD[6].lo[q] = lo; RD[6].cnt[q] = 1; RD[6].m0[q] = mm[0]; RD[6].m1[q] = mm[1];
        double a = fabs((double)lo) / 309485009821345068724781056.0; if (a > maxpar) maxpar = a;
    }
    fclose(f);
    if (np != RD[6].npairs) die("parent does not cover all legal pairs at radius 6");
    printf("parent: %u values, all legal pairs at radius 6 covered, max |x_P| %.6g\n", np, maxpar);
    /* edges: scatter */
    f = fopen(fedg, "rb"); if (!f) die("edges file");
    uint64_t ne; if (fread(&ne, 8, 1, f) != 1) die("edges header");
    uint64_t used = 0, adds = 0; uint64_t per_stage_adds[5] = {0};
    double maxterm = 0;
    for (uint64_t i = 0; i < ne; i++) {
        uint8_t hd[12]; i128 term; uint64_t mm[2];
        if (fread(hd, 12, 1, f) != 1 || fread(&term, 16, 1, f) != 1 || fread(mm, 16, 1, f) != 1) die("edge read");
        int rs = hd[0], lev = hd[1], len = hd[2], h = hd[3], mv = hd[4], st = hd[5]; uint32_t cd; memcpy(&cd, hd + 8, 4);
        if (st > SMAX) continue;
        if (rs != STAGE_R[st]) die("edge radius");
        used++;
        double a = fabs((double)term) / 309485009821345068724781056.0; if (a > maxterm) maxterm = a;
        int eL = rs - h, eR = rs - (len - 1 - h);
        if (eL < 0 || eR < 0) die("edge source does not fit its stage radius");
        Rad *R = &RD[rs]; int n = R->n;
        /* tape heights: end height */
        int hend = lev; for (int j = len - 1; j >= 0; j--) hend += ((cd >> j) & 1) ? 1 : -1;
        for (uint32_t lb = 0; lb < (1u << eL); lb++) {
            /* left letters lb (eL letters, MSB first) end at height lev: start level = lev - sum */
            int s = 2 * __builtin_popcount(lb) - eL; int wl = lev - s;
            if (!in_strip(wl, lb, eL)) continue;
            for (uint32_t rb = 0; rb < (1u << eR); rb++) {
                if (!in_strip(hend, rb, eR)) continue;
                uint64_t wc = ((uint64_t)lb << (len + eR)) | ((uint64_t)cd << eR) | rb;
                int64_t q = pair_index(R, wl, wc, mv);
                if (q < 0) die("completion window/move missing");
                if (!R->touched) { R->touched = calloc(R->npairs, 1); if (!R->touched) die("alloc touched"); }
                R->touched[q] |= (uint8_t)(1u << st);
                R->lo[q] += term; R->cnt[q] += 1;
                R->m0[q] = madd(R->m0[q], mm[0], 0); R->m1[q] = madd(R->m1[q], mm[1], 1);
                adds++; per_stage_adds[st]++;
            }
        }
    }
    fclose(f);
    printf("edges used %llu of %llu, completions added %llu (per stage %llu %llu %llu %llu %llu), max |edge term| %.6g, %.1fs\n",
           (unsigned long long)used, (unsigned long long)ne, (unsigned long long)adds,
           (unsigned long long)per_stage_adds[0], (unsigned long long)per_stage_adds[1], (unsigned long long)per_stage_adds[2],
           (unsigned long long)per_stage_adds[3], (unsigned long long)per_stage_adds[4], maxterm, now() - t0);
    fflush(stdout);
    for (int st = 0; st <= SMAX; st++) {
        Rad *R = &RD[STAGE_R[st]]; uint64_t cp = 0, cw = 0;
        if (R->touched) for (uint32_t wi = 0; wi < R->nwin; wi++) {
            int k = __builtin_popcount(R->wmask[wi]), any = 0;
            for (int j = 0; j < k; j++) if ((R->touched[R->wbase[wi] + j] >> st) & 1) { cp++; any = 1; }
            cw += any;
        }
        printf("stage %d at radius %d: changed pairs %llu, changed windows %llu\n", st, STAGE_R[st], (unsigned long long)cp, (unsigned long long)cw);
    }
    /* overflow guard: every partial sum is bounded by (number of terms) * max|term| */
    Th th[64]; pthread_t pt[64];
    memset(th, 0, sizeof th);
    for (int i = 0; i < NTH; i++) { th[i].tid = i; pthread_create(&pt[i], NULL, pass1, &th[i]); }
    for (int i = 0; i < NTH; i++) pthread_join(pt[i], NULL);
    uint64_t nonpos = 0, pairs = 0, wins = 0; double minx = 1e300; uint32_t maxcnt = 0; u128 H = 0, Mhi = 0, Mlo = 0;
    uint64_t mass[2] = {0, 0}, da[2] = {0, 0}, db[2] = {0, 0};
    for (int i = 0; i < NTH; i++) {
        nonpos += th[i].nonpos; pairs += th[i].pairs; wins += th[i].wins; if (th[i].minx < minx) minx = th[i].minx;
        if (th[i].maxcnt > maxcnt) maxcnt = th[i].maxcnt; H += th[i].H; Mhi += th[i].Mhi; Mlo += th[i].Mlo;
        for (int k = 0; k < 2; k++) { mass[k] = madd(mass[k], th[i].mass[k], k); da[k] = madd(da[k], th[i].da[k], k); db[k] = madd(db[k], th[i].db[k], k); }
    }
    double bound = (double)maxcnt * (maxterm > maxpar ? maxterm : maxpar);
    printf("pass1 %.1fs: windows %llu pairs %llu nonpositive %llu min density (approx) %.18g max terms per pair %u\n",
           now() - t0, (unsigned long long)wins, (unsigned long long)pairs, (unsigned long long)nonpos, minx, maxcnt);
    printf("overflow guard: max terms * max|term| = %.4g < 2^39 = %.4g : %s\n", bound, ldexp(1.0, 39), bound < ldexp(1.0, 39) ? "OK" : "FAIL");
    print_u128("H_lo_q120", H); print_u128("M_hi_q120", Mhi); print_u128("M_lo_q120", Mlo);
    printf("mass_mod_montgomery %llu %llu\n", (unsigned long long)mass[0], (unsigned long long)mass[1]);
    printf("drift_mod a %llu %llu b %llu %llu : %s\n", (unsigned long long)da[0], (unsigned long long)da[1],
           (unsigned long long)db[0], (unsigned long long)db[1], (da[0] | da[1] | db[0] | db[1]) ? "NONZERO" : "zero");
    fflush(stdout);
    if (ROWS) {
        memset(th, 0, sizeof th);
        for (int i = 0; i < NTH; i++) { th[i].tid = i; th[i].firstbad_key = ~0ULL; pthread_create(&pt[i], NULL, pass2, &th[i]); }
        for (int i = 0; i < NTH; i++) pthread_join(pt[i], NULL);
        uint64_t rows = 0, b0 = 0, b1 = 0, an = 0, fb = ~0ULL;
        for (int i = 0; i < NTH; i++) { rows += th[i].rows; b0 += th[i].badrows[0]; b1 += th[i].badrows[1]; an += th[i].anomalies; if (th[i].firstbad_key < fb) fb = th[i].firstbad_key; }
        printf("pass2 %.1fs: balance rows %llu, failing mod p1 %llu, failing mod p2 %llu, anomalies %llu", now() - t0,
               (unsigned long long)rows, (unsigned long long)b0, (unsigned long long)b1, (unsigned long long)an);
        if (fb != ~0ULL) printf(", first failing row key %llu (level %llu)", (unsigned long long)fb, (unsigned long long)(fb >> (RD[EVR].n + 2)));
        printf("\n");
    }
    printf("done %.1fs\n", now() - t0);
    return 0;
}
#endif
