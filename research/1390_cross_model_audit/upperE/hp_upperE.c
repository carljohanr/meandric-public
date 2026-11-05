/* hp_upperE.c -- independent checker (1390, prefix HP) for Paper E, Lemma 5 (sub-invariant weight)
 * and Lemma 6 (truncation), written from the paper's statement and the BICCERT1 format description
 * in 602 PROOF.md section 6 only.  No project checker code was read.
 *
 * Statement checked: for every state x = (U, L, c) whose stacks have at most r+1 symbols,
 *     sum over successors y of v(y) / v(x)  <=  R = P/Q,
 * with v(U,L,c) = a^{|U|+|L|} p^{E(U)+E(L)} t^{T(U)+T(L)} u(view(U), view(L), c).
 * The moves are applied literally; E and T are recounted on the whole stacks before and after;
 * views are read off the stacks.  Comparison in exact unsigned 128-bit integers after clearing
 * denominators, with a start-up bound that excludes overflow.
 *
 * Usage: hp_upperE cert.bin [threads] [dP]   (dP: subtract dP from P, for the sharpness test)
 *        hp_upperE cert.bin T trunc N          (random test of the truncation lemma on N long states)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

typedef unsigned __int128 u128;

static int r;
static uint64_t an, ad, pn, pd, tn, td, P, Q, cnt;
static uint32_t *tab;

/* A stack is (len, bits) with the TOP symbol in bit 0 and depth i in bit i. */
static inline int E_of(int len, uint64_t b) {
    int e = 0;
    for (int i = 0; i + 1 < len; i++) e += (((b >> i) & 1) == ((b >> (i + 1)) & 1));
    return e;
}
static inline int T_of(int len, uint64_t b) {
    int t = 0;
    for (int i = 0; i + 2 < len; i++) {
        int x = (b >> i) & 1, y = (b >> (i + 1)) & 1, z = (b >> (i + 2)) & 1;
        t += (x == y && y == z);
    }
    return t;
}
/* view code: a leading 1 followed by the displayed parities, top = lowest bit.
 * Shorter than r: whole stack plus bottom mark -> (1<<len) | bits.  Otherwise top r symbols -> (1<<r) | (bits & mask). */
static inline uint32_t view_code(int len, uint64_t b) {
    if (len < r) return (uint32_t)((1ull << len) | b);
    return (uint32_t)((1ull << r) | (b & ((1ull << r) - 1)));
}
static inline uint32_t swap_code(uint32_t code) {
    int len = 63 - __builtin_clzll((uint64_t)code);
    uint32_t m = (1u << len) - 1;
    return (1u << len) | ((~code) & m);
}
static inline uint64_t u_of(uint32_t cu, uint32_t cl, int c) {
    if (c) { cu = swap_code(cu); cl = swap_code(cl); }
    return tab[((uint64_t)cu << (r + 1)) | cl];
}

static u128 Aterm[2], Bterm[3], Cterm[3]; /* indexed by s in {-1,+1} -> 0,1 ; e in {-1,0,1} -> 0,1,2 */
static u128 RHSmul;
static uint64_t Pcheck;

typedef struct { int tid, nth; uint64_t nstates, nfail, nfail_hi; long double maxq; int arg_lu, arg_ll, arg_c; uint64_t arg_bu, arg_bl; long double maxf; int f_lu, f_ll, f_c; uint64_t f_bu, f_bl; } job_t;

/* evaluate one state; returns lhs (scaled) and rhs-unit (scaled without P) */
static inline void eval_state(int lu, uint64_t bu, int ll, uint64_t bl, int c, u128 *lhs, uint64_t *ux, long double *q) {
    uint64_t uxv = u_of(view_code(lu, bu), view_code(ll, bl), c);
    int EU = E_of(lu, bu), TU = T_of(lu, bu), EL = E_of(ll, bl), TL = T_of(ll, bl);
    u128 s = 0; long double fq = 0;
    long double fa = (long double)an / ad, fp = (long double)pn / pd, ft = (long double)tn / td;
    /* push on U */
    {
        int nl = lu + 1; uint64_t nb = (bu << 1) | (uint64_t)c;
        int dE = E_of(nl, nb) - EU, dT = T_of(nl, nb) - TU;
        uint64_t uy = u_of(view_code(nl, nb), view_code(ll, bl), 1 - c);
        s += Aterm[1] * Bterm[dE + 1] * Cterm[dT + 1] * (u128)Q * uy;
        fq += fa * (dE ? fp : 1) * (dT ? ft : 1) * (long double)uy / uxv;
    }
    /* push on L */
    {
        int nl = ll + 1; uint64_t nb = (bl << 1) | (uint64_t)c;
        int dE = E_of(nl, nb) - EL, dT = T_of(nl, nb) - TL;
        uint64_t uy = u_of(view_code(lu, bu), view_code(nl, nb), 1 - c);
        s += Aterm[1] * Bterm[dE + 1] * Cterm[dT + 1] * (u128)Q * uy;
        fq += fa * (dE ? fp : 1) * (dT ? ft : 1) * (long double)uy / uxv;
    }
    /* pop U: allowed if nonempty and top == 1-c */
    if (lu >= 1 && (int)(bu & 1) == 1 - c) {
        int nl = lu - 1; uint64_t nb = bu >> 1;
        int dE = E_of(nl, nb) - EU, dT = T_of(nl, nb) - TU;
        uint64_t uy = u_of(view_code(nl, nb), view_code(ll, bl), 1 - c);
        s += Aterm[0] * Bterm[dE + 1] * Cterm[dT + 1] * (u128)Q * uy;
        fq += (1 / fa) * (dE ? 1 / fp : 1) * (dT ? 1 / ft : 1) * (long double)uy / uxv;
    }
    if (ll >= 1 && (int)(bl & 1) == 1 - c) {
        int nl = ll - 1; uint64_t nb = bl >> 1;
        int dE = E_of(nl, nb) - EL, dT = T_of(nl, nb) - TL;
        uint64_t uy = u_of(view_code(lu, bu), view_code(nl, nb), 1 - c);
        s += Aterm[0] * Bterm[dE + 1] * Cterm[dT + 1] * (u128)Q * uy;
        fq += (1 / fa) * (dE ? 1 / fp : 1) * (dT ? 1 / ft : 1) * (long double)uy / uxv;
    }
    *lhs = s; *ux = uxv; *q = fq;
}

static void *worker(void *arg) {
    job_t *j = (job_t *)arg;
    int K = r + 1;
    uint64_t idx = 0;
    for (int lu = 0; lu <= K; lu++)
        for (uint64_t bu = 0; bu < (1ull << lu); bu++, idx++) {
            if ((int)(idx % j->nth) != j->tid) continue;
            for (int ll = 0; ll <= K; ll++)
                for (uint64_t bl = 0; bl < (1ull << ll); bl++)
                    for (int c = 0; c < 2; c++) {
                        u128 lhs; uint64_t ux; long double q;
                        eval_state(lu, bu, ll, bl, c, &lhs, &ux, &q);
                        u128 rhs = (u128)Pcheck * RHSmul * ux;
                        u128 rhs_hi = (u128)P * RHSmul * ux;
                        j->nstates++;
                        if (lhs > rhs) j->nfail++;
                        if (lhs > rhs_hi) j->nfail_hi++;
                        if (q > j->maxq) { j->maxq = q; j->arg_lu = lu; j->arg_bu = bu; j->arg_ll = ll; j->arg_bl = bl; j->arg_c = c; }
                        if (lu >= r && ll >= r && q > j->maxf) { j->maxf = q; j->f_lu = lu; j->f_bu = bu; j->f_ll = ll; j->f_bl = bl; j->f_c = c; }
                    }
        }
    return NULL;
}

static void print_stack(int len, uint64_t b) { /* bottom first */
    if (len == 0) { printf("(empty)"); return; }
    for (int i = len - 1; i >= 0; i--) putchar('0' + (int)((b >> i) & 1));
}

static uint64_t rng = 88172645463325252ull;
static inline uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 2; }
    int nth = argc > 2 ? atoi(argv[2]) : 4;
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror("open"); return 2; }
    char magic[8];
    if (fread(magic, 1, 8, f) != 8 || memcmp(magic, "BICCERT1", 8)) { printf("REJECT: bad magic\n"); return 1; }
    int32_t r32; if (fread(&r32, 4, 1, f) != 1) { printf("REJECT: short header\n"); return 1; }
    r = r32;
    if (r < 2 || r > 14) { printf("REJECT: radius %d out of range\n", r); return 1; }
    uint64_t h[9];
    if (fread(h, 8, 9, f) != 9) { printf("REJECT: short header\n"); return 1; }
    an = h[0]; ad = h[1]; pn = h[2]; pd = h[3]; tn = h[4]; td = h[5]; P = h[6]; Q = h[7]; cnt = h[8];
    printf("r=%d a=%llu/%llu p=%llu/%llu t=%llu/%llu R=%llu/%llu count=%llu\n", r,
           (unsigned long long)an, (unsigned long long)ad, (unsigned long long)pn, (unsigned long long)pd,
           (unsigned long long)tn, (unsigned long long)td, (unsigned long long)P, (unsigned long long)Q, (unsigned long long)cnt);
    if (!an || !ad || !pn || !pd || !tn || !td || !P || !Q) { printf("REJECT: zero parameter\n"); return 1; }
    if (cnt != (1ull << (2 * (r + 1)))) { printf("REJECT: count != 4^(r+1)\n"); return 1; }
    tab = malloc(cnt * 4);
    if (fread(tab, 4, cnt, f) != cnt) { printf("REJECT: truncated table\n"); return 1; }
    int extra = fgetc(f);
    if (extra != EOF) { printf("REJECT: trailing bytes\n"); return 1; }
    fclose(f);
    uint64_t nV = (1ull << (r + 1)) - 1, umax = 0;
    for (uint64_t x = 1; x <= nV; x++)
        for (uint64_t y = 1; y <= nV; y++) {
            uint32_t e = tab[(x << (r + 1)) | y];
            if (e == 0) { printf("REJECT: zero entry at used position (%llu,%llu)\n", (unsigned long long)x, (unsigned long long)y); return 1; }
            if (e > umax) umax = e;
        }
    printf("used entries: %llu, all positive, max %llu\n", (unsigned long long)(nV * nV), (unsigned long long)umax);
    Aterm[1] = (u128)an * an; Aterm[0] = (u128)ad * ad;
    Bterm[2] = (u128)pn * pn; Bterm[1] = (u128)pn * pd; Bterm[0] = (u128)pd * pd;
    Cterm[2] = (u128)tn * tn; Cterm[1] = (u128)tn * td; Cterm[0] = (u128)td * td;
    RHSmul = (u128)an * ad * pn * pd * tn * td;
    /* overflow guard: bound every quantity by 2^127 using doubles with a generous margin */
    long double mA = (long double)(an > ad ? an : ad), mB = (long double)(pn > pd ? pn : pd), mC = (long double)(tn > td ? tn : td);
    long double lhsmax = 4.0L * mA * mA * mB * mB * mC * mC * (long double)Q * (long double)umax;
    long double rhsmax = (long double)P * (long double)an * ad * pn * pd * tn * td * (long double)umax;
    long double lim = 1.0L; for (int i = 0; i < 126; i++) lim *= 2;  /* 2^126, margin 2x */
    printf("overflow guard: lhs<=%.3Le rhs<=%.3Le limit 2^126=%.3Le\n", lhsmax, rhsmax, lim);
    if (lhsmax >= lim || rhsmax >= lim) { printf("ABORT: overflow possible\n"); return 3; }

    if (argc > 4 && !strcmp(argv[3], "trunc")) {
        /* random test of the truncation lemma: long states vs. their truncation to r+1 symbols */
        long N = atol(argv[4]); long bad = 0;
        for (long i = 0; i < N; i++) {
            int lu = (int)(xr() % 40), ll = (int)(xr() % 40), c = (int)(xr() & 1);
            uint64_t bu = xr() & ((1ull << lu) - 1), bl = xr() & ((1ull << ll) - 1);
            u128 l1, l2; uint64_t u1, u2; long double q1, q2;
            eval_state(lu, bu, ll, bl, c, &l1, &u1, &q1);
            int K = r + 1;
            int tlu = lu < K ? lu : K, tll = ll < K ? ll : K;
            eval_state(tlu, bu & ((1ull << tlu) - 1), tll, bl & ((1ull << tll) - 1), c, &l2, &u2, &q2);
            /* l/u equal as rationals: l1*u2 == l2*u1, compare in 256-bit via splitting is avoided: u1==u2 is expected */
            if (u1 != u2 || l1 != l2) bad++;
        }
        printf("truncation test: %ld random states (stacks up to 39 symbols), %ld mismatches\n", N, bad);
        return bad ? 1 : 0;
    }

    uint64_t dP = 0;
    const char *dmg = NULL;
    if (argc > 3 && !strcmp(argv[3], "dmg")) dmg = argv[4]; else if (argc > 3) dP = strtoull(argv[3], 0, 10);
    Pcheck = P - dP;
    pthread_t th[64]; job_t jobs[64];
    uint64_t ns, nf, nfh; long double mq; int best;
    for (int pass = 0; pass < (dmg ? 2 : 1); pass++) {
    for (int i = 0; i < nth; i++) { memset(&jobs[i], 0, sizeof(job_t)); jobs[i].tid = i; jobs[i].nth = nth; pthread_create(&th[i], 0, worker, &jobs[i]); }
    ns = 0; nf = 0; nfh = 0; mq = 0; best = 0; long double mf = 0; int bestf = 0;
    for (int i = 0; i < nth; i++) {
        pthread_join(th[i], 0); ns += jobs[i].nstates; nf += jobs[i].nfail; nfh += jobs[i].nfail_hi;
        if (jobs[i].maxq > mq) { mq = jobs[i].maxq; best = i; }
        if (jobs[i].maxf > mf) { mf = jobs[i].maxf; bestf = i; }
    }
    if (pass == 0) {
        printf("tightest state with two full views: q=%.15Lf U=", mf); print_stack(jobs[bestf].f_lu, jobs[bestf].f_bu);
        printf(" L="); print_stack(jobs[bestf].f_ll, jobs[bestf].f_bl); printf(" c=%d\n", jobs[bestf].f_c);
    }
    if (dmg && pass == 0) {
        job_t *J = &jobs[best]; int lu = J->arg_lu, ll = J->arg_ll, c = J->arg_c; uint64_t bu = J->arg_bu, bl = J->arg_bl;
        if (!strcmp(dmg, "fullview")) { J = &jobs[bestf]; lu = J->f_lu; ll = J->f_ll; c = J->f_c; bu = J->f_bu; bl = J->f_bl; }
        uint32_t cu = view_code(lu, bu), cl = view_code(ll, bl);
        if (!strcmp(dmg, "succ_up")) { /* successor by push on L */ cl = view_code(ll + 1, (bl << 1) | (uint64_t)c); c = 1 - c; }
        if (c) { cu = swap_code(cu); cl = swap_code(cl); }
        uint64_t idx = ((uint64_t)cu << (r + 1)) | cl;
        uint32_t old = tab[idx];
        if (!strcmp(dmg, "tight_down") || !strcmp(dmg, "fullview")) tab[idx] = (uint32_t)(old - old / 1000);
        else if (!strcmp(dmg, "succ_up")) tab[idx] = (uint32_t)((uint64_t)old + old / 1000 > 0xffffffffull ? 0xffffffffull : (uint64_t)old + old / 1000);
        else if (!strcmp(dmg, "zero")) tab[idx] = 0;
        printf("DAMAGE %s: entry %llu changed %u -> %u\n", dmg, (unsigned long long)idx, old, tab[idx]);
        if (tab[idx] == 0) { printf("REJECT: zero entry at used position\nRESULT: REJECT\n"); return 1; }
        if (!strcmp(dmg, "succ_up") && tab[idx] == old) { printf("damage saturated\n"); return 4; }
    }
    }
    uint64_t expect = 2 * ((1ull << (r + 2)) - 1) * ((1ull << (r + 2)) - 1);
    printf("states checked: %llu (expected 2(2^(r+2)-1)^2 = %llu)\n", (unsigned long long)ns, (unsigned long long)expect);
    printf("violations at R=(P-%llu)/Q: %llu ; at stated R: %llu\n", (unsigned long long)dP, (unsigned long long)nf, (unsigned long long)nfh);
    printf("largest quotient (long double): %.15Lf at U=", mq);
    print_stack(jobs[best].arg_lu, jobs[best].arg_bu); printf(" L="); print_stack(jobs[best].arg_ll, jobs[best].arg_bl);
    printf(" c=%d\n", jobs[best].arg_c);
    if (ns != expect) { printf("RESULT: ERROR state count\n"); return 1; }
    printf("RESULT: %s\n", nf == 0 ? "ACCEPT" : "REJECT");
    return nf == 0 ? 0 : 1;
}
