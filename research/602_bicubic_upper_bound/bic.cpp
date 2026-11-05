// 602: window-weight solver for the two-stack transfer of Hamiltonian cycles on bicubic planar maps.
//
// State (U, L, c): two stacks of parities, c = parity of the current position. Moves at parity c:
//   push c on U, push c on L, pop U if its top is 1-c, pop L if its top is 1-c; then c flips.
// Weight: v = a^(|U|+|L|) p^(E(U)+E(L)) t^(T(U)+T(L)) u(view_r(U), view_r(L), c),
//   E = adjacent equal pairs in a stack, T = adjacent equal triples.
// A view is coded as an integer: a leading 1 followed by the displayed parities (top = lowest bit). Codes below
// 2^r display the whole stack (fewer than r parities above the bottom mark); codes 2^r .. 2^(r+1)-1 display r
// parities and nothing below them.
// The table is stored for c = 0 only: u(x, y, 1) = u(comp x, comp y, 0), comp = exchange of the two parities.
// The solver also keeps the table symmetric under exchange of the two stacks. These are choices of the search;
// any positive table gives a bound (see PROOF.md); the verifiers check the expanded table.
//
// modes:
//   bic solve r a p t [in=F] [out=F] [tol=1e-10] [maxit=N] [maxsec=S] [threads=2] [float=1] [hb=0.75] [quiet=1]
//   bic cert  r an ad pn pd tn td Q in=F out=F [float=1] [threads=2] [scale=4294967295]
//   bic check certfile [threads=2]        exact check of every row (the solver's own row formulas)
// Requires r >= 3 (so that an equal triple at the top never involves the revealed symbol).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <algorithm>
#include <map>
typedef uint64_t u64; typedef uint32_t u32; typedef unsigned __int128 u128;

static int r_; static u64 M_, FULL_;                    // M = 2^(r+1), FULL = 2^r
static inline int blen(u64 code) { return 63 - __builtin_clzll(code); }      // number of displayed parities
static inline u64 comp(u64 code) { return code ^ ((1ull << blen(code)) - 1); }
static inline u64 pushc(u64 code, int c) { return code < FULL_ ? ((code << 1) | (u64)c) : (FULL_ | (((code << 1) | (u64)c) & (FULL_ - 1))); }

// Per view and parity c: everything a move on that stack needs.
struct ViewMove {
    u64 pushDest; int pushE, pushT;                     // push c
    int popN; u64 popDest[3]; int popE, popT;           // pop (popN = 0: not available)
};
static ViewMove viewMove(u64 code, int c) {
    ViewMove m; int n = blen(code);
    int top = n >= 1 ? (int)(code & 1) : -1, sec = n >= 2 ? (int)((code >> 1) & 1) : -1, thi = n >= 3 ? (int)((code >> 2) & 1) : -1;
    m.pushDest = pushc(code, c); m.pushE = (top == c); m.pushT = (top == c && sec == c);
    m.popN = 0; m.popE = m.popT = 0; m.popDest[0] = m.popDest[1] = m.popDest[2] = 0;
    if (top >= 0 && top != c) {
        m.popE = (sec == top); m.popT = (sec == top && thi == top);
        if (code < FULL_) { m.popN = 1; m.popDest[0] = code >> 1; }
        else { u64 b = (code & (FULL_ - 1)) >> 1; m.popN = 3; m.popDest[0] = code >> 1; m.popDest[1] = FULL_ | b; m.popDest[2] = FULL_ | (FULL_ >> 1) | b; }
    }
    return m;
}

static std::map<std::string, std::string> kv;
static double argd(const char* k, double d) { auto it = kv.find(k); return it == kv.end() ? d : atof(it->second.c_str()); }
static std::string args(const char* k) { auto it = kv.find(k); return it == kv.end() ? std::string() : it->second; }
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

template <class FT> struct Solver {
    std::vector<FT> u, Tu; int threads;
    // column data (stack L), per code y, for c = 0 at the source and the destination read at c = 1 through comp
    std::vector<u32> cY, pushY, popY0, popY1, popY2; std::vector<double> wPushY, wPopY; std::vector<unsigned char> popNY;
    double a, p, t;
    void setup(double a_, double p_, double t_) {
        a = a_; p = p_; t = t_;
        cY.assign(M_, 0); pushY.assign(M_, 0); popY0.assign(M_, 0); popY1.assign(M_, 0); popY2.assign(M_, 0);
        wPushY.assign(M_, 0); wPopY.assign(M_, 0); popNY.assign(M_, 0);
        for (u64 y = 1; y < M_; y++) {
            ViewMove m = viewMove(y, 0);
            cY[y] = (u32)comp(y); pushY[y] = (u32)comp(m.pushDest); wPushY[y] = a * (m.pushE ? p : 1) * (m.pushT ? t : 1);
            popNY[y] = (unsigned char)m.popN;
            if (m.popN) { wPopY[y] = (1 / a) / (m.popE ? p : 1) / (m.popT ? t : 1); popY0[y] = (u32)comp(m.popDest[0]); popY1[y] = popY0[y]; popY2[y] = popY0[y];
                if (m.popN == 3) { popY1[y] = (u32)comp(m.popDest[1]); popY2[y] = (u32)comp(m.popDest[2]); } }
        }
    }
    // One step on the upper triangle (y >= x), mirrored. Reads u (= u_k); the array Tu holds u_(k-1) on entry and
    // u_(k+1) = (T u_k / Rn) * (u_k / u_(k-1))^beta on exit (beta = 0: plain power step). hi, lo = extremes of T u_k / u_k.
    // Any positive update is legitimate: the bound comes from the exact check of the final integer table.
    void apply(double& hi, double& lo, double beta, double Rn, double theta = 1.0) {
        std::atomic<u64> next(1); std::vector<double> his(threads, 0), los(threads, 1e300);
        auto work = [&](int tid) {
            double h = 0, l = 1e300; const FT* S = u.data(); double inv = 1 / Rn;
            for (;;) {
                u64 x = next.fetch_add(1); if (x >= M_) break;
                ViewMove m = viewMove(x, 0);
                const FT* rowPush = S + comp(m.pushDest) * M_; double wPush = wPushY[x];
                const FT* rowSame = S + comp(x) * M_;
                const FT* rp0 = S, *rp1 = S, *rp2 = S; double wPop = 0;
                if (m.popN) { wPop = wPopY[x]; rp0 = S + comp(m.popDest[0]) * M_; rp1 = rp2 = rp0;
                    if (m.popN == 3) { rp1 = S + comp(m.popDest[1]) * M_; rp2 = S + comp(m.popDest[2]) * M_; } }
                for (u64 y = x; y < M_; y++) {
                    u32 cy = cY[y];
                    double v = wPush * rowPush[cy] + wPushY[y] * rowSame[pushY[y]];
                    if (m.popN) { double z = rp0[cy]; if (m.popN == 3) { z = std::max(z, (double)rp1[cy]); z = std::max(z, (double)rp2[cy]); } v += wPop * z; }
                    if (popNY[y]) { double z = rowSame[popY0[y]]; if (popNY[y] == 3) { z = std::max(z, (double)rowSame[popY1[y]]); z = std::max(z, (double)rowSame[popY2[y]]); } v += wPopY[y] * z; }
                    double cur = S[x * M_ + y], rat = v / cur; if (rat > h) h = rat; if (rat < l) l = rat;
                    double nv = v * inv;
                    if (theta != 1.0) nv = cur * std::pow(nv / cur, theta);
                    if (beta > 0) { double q = cur / (double)Tu[x * M_ + y]; nv *= (beta == 0.75 ? std::sqrt(q) * std::sqrt(std::sqrt(q)) : beta == 0.5 ? std::sqrt(q) : std::pow(q, beta)); }
                    Tu[x * M_ + y] = (FT)nv; Tu[y * M_ + x] = (FT)nv;
                }
            }
            his[tid] = h; los[tid] = l;
        };
        std::vector<std::thread> th; for (int i = 0; i < threads; i++) th.emplace_back(work, i); for (auto& q : th) q.join();
        hi = *std::max_element(his.begin(), his.end()); lo = *std::min_element(los.begin(), los.end());
    }
};

// table files: magic "BICTAB1\0", int32 r, int32 bytes per entry (4 float, 8 double), then M*M entries (c = 0)
template <class FT> static bool loadTable(const std::string& f, std::vector<FT>& u) {
    FILE* fp = fopen(f.c_str(), "rb"); if (!fp) return false; char mg[8]; int rr, by;
    if (fread(mg, 1, 8, fp) != 8 || memcmp(mg, "BICTAB1", 7) || fread(&rr, 4, 1, fp) != 1 || fread(&by, 4, 1, fp) != 1) { fclose(fp); return false; }
    u64 Ms = 1ull << (rr + 1), Fs = 1ull << rr; std::vector<double> s(Ms * Ms);
    if (by == 8) { if (fread(s.data(), 8, Ms * Ms, fp) != Ms * Ms) { fclose(fp); return false; } }
    else { std::vector<float> g(Ms * Ms); if (fread(g.data(), 4, Ms * Ms, fp) != Ms * Ms) { fclose(fp); return false; } for (u64 i = 0; i < Ms * Ms; i++) s[i] = g[i]; }
    fclose(fp);
    if (rr == r_) { for (u64 i = 0; i < M_ * M_; i++) u[i] = (FT)s[i]; return true; }
    if (rr > r_) return false;
    // lift from a smaller radius: forget the deepest displayed parities
    auto trunc = [&](u64 code) { while (code >= 2 * Fs) { int n = blen(code); code = (1ull << (n - 1)) | (code & ((1ull << (n - 1)) - 1)); } return code; };
    std::vector<u64> tr(M_); for (u64 x = 1; x < M_; x++) tr[x] = trunc(x);
    for (u64 x = 1; x < M_; x++) for (u64 y = 1; y < M_; y++) u[x * M_ + y] = (FT)s[tr[x] * Ms + tr[y]];
    return true;
}
template <class FT> static void saveTable(const std::string& f, const std::vector<FT>& u) {
    FILE* fp = fopen(f.c_str(), "wb"); int by = sizeof(FT); fwrite("BICTAB1\0", 1, 8, fp); fwrite(&r_, 4, 1, fp); fwrite(&by, 4, 1, fp); fwrite(u.data(), by, M_ * M_, fp); fclose(fp);
}

template <class FT> static int runSolve(double a, double p, double t) {
    Solver<FT> S; S.threads = (int)argd("threads", 2); S.u.assign(M_ * M_, 0); S.Tu.assign(M_ * M_, 0);
    for (u64 x = 1; x < M_; x++) for (u64 y = 1; y < M_; y++) S.u[x * M_ + y] = 1;
    if (!args("in").empty() && !loadTable(args("in"), S.u)) { fprintf(stderr, "cannot read %s\n", args("in").c_str()); return 1; }
    S.setup(a, p, t);
    double tol = argd("tol", 1e-10), maxsec = argd("maxsec", 1100), hb = argd("hb", 0.75), theta = argd("theta", 1.0); int maxit = (int)argd("maxit", 100000), quiet = (int)argd("quiet", 0);
    double t0 = now(), hi = 0, lo = 0, Rn = 3.2, bestSpread = 1e300; int it = 0, conv = 0, cool = 1, resets = 0;
    for (it = 1; it <= maxit; it++) {
        // momentum (heavy ball in the logarithm) is switched off for the first step and for 5 steps after the spread
        // has grown to 4 times its best value
        S.apply(hi, lo, cool > 0 ? 0.0 : hb, Rn, theta); if (cool > 0) cool--;
        std::swap(S.u, S.Tu);                           // now Tu = u_k (the table hi, lo refer to), u = u_(k+1)
        if (!(hi >= lo) || !std::isfinite(hi) || !(lo > 0)) { printf("DIVERGED at iteration %d\n", it); return 1; }
        double spread = hi - lo; Rn = std::sqrt(hi * lo);
        if (!quiet && (it % 10 == 0 || it <= 3)) { printf("   it %d: hi %.10f lo %.10f spread %.3e (%.1f s)\n", it, hi, lo, spread, now() - t0); fflush(stdout); }
        if (spread < tol * hi) { conv = 1; break; }
        if (now() - t0 > maxsec) break;
        if (spread < bestSpread) bestSpread = spread; else if (spread > 4 * bestSpread && cool == 0 && it > 20) { cool = 5; resets++; bestSpread = spread; }
    }
    // save u_k: its ratios are within [lo, hi]
    if (!args("out").empty()) saveTable(args("out"), S.Tu);
    if (resets) printf("   momentum resets: %d\n", resets);
    printf("r=%d a=%.6f p=%.6f t=%.6f: R_float %.9f (hi %.9f lo %.9f) R^2 %.6f iterations %d converged %d entries %llu, %.1f s\n",
        r_, a, p, t, hi, hi, lo, hi * hi, it, conv, (unsigned long long)((M_ - 1) * (M_ - 1)), now() - t0);
    fflush(stdout); return 0;
}

// ---- exact certificate ----
// file: magic "BICCERT1", int32 r, uint64 an, ad, pn, pd, tn, td, P, Q, uint64 count = (2^(r+1))^2, then count uint32:
// entry [x * 2^(r+1) + y] = u(x, y, 0) for view codes x, y >= 1 (entries with x = 0 or y = 0 are 0 and unused).
struct Cert { int r; u64 an, ad, pn, pd, tn, td, P, Q; std::vector<u32> tab; };
static bool readCert(const std::string& f, Cert& c) {
    FILE* fp = fopen(f.c_str(), "rb"); if (!fp) return false; char mg[8]; u64 h[9];
    if (fread(mg, 1, 8, fp) != 8 || memcmp(mg, "BICCERT1", 8) || fread(&c.r, 4, 1, fp) != 1 || fread(h, 8, 9, fp) != 9) { fclose(fp); return false; }
    c.an = h[0]; c.ad = h[1]; c.pn = h[2]; c.pd = h[3]; c.tn = h[4]; c.td = h[5]; c.P = h[6]; c.Q = h[7];
    u64 Ms = 1ull << (c.r + 1); if (h[8] != Ms * Ms) { fclose(fp); return false; }
    c.tab.resize(Ms * Ms); bool ok = fread(c.tab.data(), 4, Ms * Ms, fp) == Ms * Ms; char x; ok = ok && fread(&x, 1, 1, fp) == 0; fclose(fp); return ok;
}
static void writeCert(const std::string& f, const Cert& c) {
    FILE* fp = fopen(f.c_str(), "wb"); fwrite("BICCERT1", 1, 8, fp); fwrite(&c.r, 4, 1, fp);
    u64 h[9] = { c.an, c.ad, c.pn, c.pd, c.tn, c.td, c.P, c.Q, (u64)c.tab.size() }; fwrite(h, 8, 9, fp); fwrite(c.tab.data(), 4, c.tab.size(), fp); fclose(fp);
}
// Scaled row inequality. Multiply  a p^e t^f u + ... + (1/a) p^-e t^-f u + ... <= (P/Q) u0  by an ad pn pd tn td Q:
//   push term: an^2 * Ep(e) * Et(f) * u,   Ep(0) = pn pd, Ep(1) = pn^2;  Et likewise
//   pop  term: ad^2 * Fp(e) * Ft(f) * u,   Fp(0) = pn pd, Fp(1) = pd^2;  Ft likewise
//   right side: P * an ad pn pd tn td * u0;  left side is multiplied by Q.
struct Exact {
    const Cert& c; u64 an2, ad2, Ep[2], Et[2], Fp[2], Ft[2]; u128 rhsK;
    Exact(const Cert& c_) : c(c_) {
        an2 = c.an * c.an; ad2 = c.ad * c.ad; Ep[0] = Fp[0] = c.pn * c.pd; Ep[1] = c.pn * c.pn; Fp[1] = c.pd * c.pd;
        Et[0] = Ft[0] = c.tn * c.td; Et[1] = c.tn * c.tn; Ft[1] = c.td * c.td;
        rhsK = (u128)c.an * c.ad * c.pn * c.pd * c.tn * c.td;
    }
    bool sane() const {                                  // excludes overflow of the 128-bit arithmetic used below
        u64 lim = 100000; if (c.an > lim || c.ad > lim || c.pn > lim || c.pd > lim || c.tn > lim || c.td > lim) return false;
        if (!c.an || !c.ad || !c.pn || !c.pd || !c.tn || !c.td || !c.P || !c.Q) return false;
        long double coef = (long double)std::max(an2, ad2) * (long double)std::max(Ep[1], Fp[1]) * (long double)std::max(Et[1], Ft[1]);
        coef = std::max(coef, (long double)std::max(an2, ad2) * (long double)Ep[0] * (long double)Et[0]);
        long double lhs = 4 * coef * 4294967296.0L * (long double)c.Q;          // four terms, entries < 2^32, then times Q
        long double rhs = (long double)rhsK * 4294967296.0L * (long double)c.P;
        long double cap = 1.0e38L;                                              // 2^128 = 3.4e38
        return coef < 1.8e19L && lhs < cap && rhs < cap;
    }
};
static int runCheck(const Cert& c, int threads, bool findP, u64* minP) {
    r_ = c.r; M_ = 1ull << (r_ + 1); FULL_ = 1ull << r_;
    Exact E(c);
    const u32* S = c.tab.data();
    for (u64 x = 1; x < M_; x++) for (u64 y = 1; y < M_; y++) if (!S[x * M_ + y]) { printf("zero entry\n"); return 1; }
    std::vector<ViewMove> vm0(M_), vm1(M_); for (u64 x = 1; x < M_; x++) { vm0[x] = viewMove(x, 0); vm1[x] = viewMove(x, 1); }
    std::vector<u32> cmp(M_); for (u64 x = 1; x < M_; x++) cmp[x] = (u32)comp(x);
    std::atomic<u64> next(1), bad(0), rows(0), digest(0); std::vector<u64> pmax(threads, 0);
    auto val = [&](u64 x, u64 y, int cc) -> u64 { return cc == 0 ? S[x * M_ + y] : S[(u64)cmp[x] * M_ + cmp[y]]; };
    const u64 FNV = 1099511628211ull;
    auto mix = [&](u64 h, u64 item) -> u64 { return (h ^ item) * FNV; };
    auto work = [&](int tid) {
        u64 nb = 0, nr = 0, pm = 0, dg = 0;
        for (;;) {
            u64 x = next.fetch_add(1); if (x >= M_) break;
            for (int cc = 0; cc < 2; cc++) {
                const ViewMove& mx = cc ? vm1[x] : vm0[x]; int nc = 1 - cc;
                for (u64 y = 1; y < M_; y++) {
                    const ViewMove& my = cc ? vm1[y] : vm0[y];
                    u128 push = (u128)(E.an2 * E.Ep[mx.pushE]) * E.Et[mx.pushT] * val(mx.pushDest, y, nc)
                              + (u128)(E.an2 * E.Ep[my.pushE]) * E.Et[my.pushT] * val(x, my.pushDest, nc);
                    u128 kx = mx.popN ? (u128)(E.ad2 * E.Fp[mx.popE]) * E.Ft[mx.popT] : 0, ky = my.popN ? (u128)(E.ad2 * E.Fp[my.popE]) * E.Ft[my.popT] : 0;
                    u128 rhs = E.rhsK * val(x, y, cc);
                    // one row per pair of revealed symbols (three choices under a view that displays r parities)
                    int nx = x >= FULL_ ? 3 : 1, ny = y >= FULL_ ? 3 : 1;
                    u64 h0 = mix(1469598103934665603ull, (x * M_ + y) * 2 + cc);
                    h0 = mix(mix(mix(mix(h0, 0), (mx.pushDest * M_ + y) * 2 + nc), mx.pushE + 1), mx.pushT + 1);
                    h0 = mix(mix(mix(mix(h0, 1), (x * M_ + my.pushDest) * 2 + nc), my.pushE + 1), my.pushT + 1);
                    for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) {
                        u128 lhs = push; u64 h = h0;
                        if (mx.popN) { lhs += kx * val(mx.popDest[i], y, nc); h = mix(mix(mix(mix(h, 2), (mx.popDest[i] * M_ + y) * 2 + nc), 1 - mx.popE), 1 - mx.popT); }
                        if (my.popN) { lhs += ky * val(x, my.popDest[j], nc); h = mix(mix(mix(mix(h, 3), (x * M_ + my.popDest[j]) * 2 + nc), 1 - my.popE), 1 - my.popT); }
                        nr++; dg += h;
                        if (findP) { u128 q = (lhs * c.Q + rhs - 1) / rhs; if ((u64)q > pm) pm = (u64)q; }
                        else if (lhs * c.Q > rhs * c.P) nb++;
                    }
                }
            }
        }
        bad += nb; rows += nr; digest += dg; pmax[tid] = pm;
    };
    std::vector<std::thread> th; for (int i = 0; i < threads; i++) th.emplace_back(work, i); for (auto& q : th) q.join();
    if (findP) *minP = *std::max_element(pmax.begin(), pmax.end());
    else printf("rows %llu, violated %llu, row digest %016llx\n", (unsigned long long)rows.load(), (unsigned long long)bad.load(), (unsigned long long)digest.load());
    return bad.load() ? 1 : 0;
}

template <class FT> static int runCert(int argc, char** argv) {
    Cert c; c.r = r_; c.an = strtoull(argv[3], 0, 10); c.ad = strtoull(argv[4], 0, 10); c.pn = strtoull(argv[5], 0, 10); c.pd = strtoull(argv[6], 0, 10);
    c.tn = strtoull(argv[7], 0, 10); c.td = strtoull(argv[8], 0, 10); c.Q = strtoull(argv[9], 0, 10); c.P = 1;
    std::vector<FT> u(M_ * M_, 0); if (!loadTable(args("in"), u)) { fprintf(stderr, "cannot read table\n"); return 1; }
    double mx = 0; for (auto z : u) if (z > mx) mx = z; double scale = argd("scale", 4294967295.0) / mx;
    c.tab.assign(M_ * M_, 0); u32 mn = 0xffffffffu;
    for (u64 x = 1; x < M_; x++) for (u64 y = 1; y < M_; y++) { double z = std::floor((double)u[x * M_ + y] * scale + 0.5); if (z < 1) z = 1; if (z > 4294967295.0) z = 4294967295.0; c.tab[x * M_ + y] = (u32)z; if ((u32)z < mn) mn = (u32)z; }
    std::vector<FT>().swap(u);
    Exact E(c); if (!E.sane()) { fprintf(stderr, "parameters too large for 128-bit arithmetic\n"); return 1; }
    u64 P = 0; if (runCheck(c, (int)argd("threads", 2), true, &P)) return 1;
    c.P = P; writeCert(args("out"), c);
    printf("certificate r=%d a=%llu/%llu p=%llu/%llu t=%llu/%llu: smallest P = %llu with Q = %llu, R = %.12f, R^2 = %.12f, smallest entry %u\n", c.r,
        (unsigned long long)c.an, (unsigned long long)c.ad, (unsigned long long)c.pn, (unsigned long long)c.pd, (unsigned long long)c.tn, (unsigned long long)c.td,
        (unsigned long long)P, (unsigned long long)c.Q, (double)P / c.Q, ((double)P / c.Q) * ((double)P / c.Q), mn);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: see the head of bic.cpp\n"); return 1; }
    std::string mode = argv[1];
    for (int i = 2; i < argc; i++) { const char* e = strchr(argv[i], '='); if (e) kv[std::string(argv[i], e - argv[i])] = e + 1; }
    if (mode == "check") {
        Cert c; if (!readCert(argv[2], c)) { printf("cannot read certificate (format error)\n"); return 1; }
        Exact E(c); if (!E.sane()) { printf("parameters out of range\n"); return 1; }
        printf("certificate r=%d a=%llu/%llu p=%llu/%llu t=%llu/%llu P/Q=%llu/%llu\n", c.r, (unsigned long long)c.an, (unsigned long long)c.ad, (unsigned long long)c.pn,
            (unsigned long long)c.pd, (unsigned long long)c.tn, (unsigned long long)c.td, (unsigned long long)c.P, (unsigned long long)c.Q);
        int rc = runCheck(c, (int)argd("threads", 2), false, nullptr); printf(rc ? "REJECTED\n" : "ACCEPTED\n"); return rc;
    }
    r_ = atoi(argv[2]); if (r_ < 3 || r_ > 14) { fprintf(stderr, "r must be 3..14\n"); return 1; }
    M_ = 1ull << (r_ + 1); FULL_ = 1ull << r_;
    bool fl = argd("float", 0) != 0;
    if (mode == "solve") { double a = atof(argv[3]), p = atof(argv[4]), t = atof(argv[5]); return fl ? runSolve<float>(a, p, t) : runSolve<double>(a, p, t); }
    if (mode == "cert") { if (argc < 10) return 1; return fl ? runCert<float>(argc, argv) : runCert<double>(argc, argv); }
    return 1;
}
