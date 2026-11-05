"""1426 (JB): read parent_r6.bin (format IWPAR001, SPEC Section 4) and check C2 and C3 exactly.

C3: the records are exactly the legal pairs at radius 6, in canonical order; sum N = D; sum N*Da = sum N*Db = 0.
C2: all 220,110 balance rows (S_6) of SPEC 7.1 hold exactly (integer arithmetic after clearing denominators).
Then writes work/parent_c.bin for the C checker: per record level(u8) move(u8) pad(u16) code(u32)
lo(int128, floor(x_P * 2^88)) modp1(u64) modp2(u64) [Montgomery form of x_P mod p].

Usage: python3 parent.py <parent_r6.bin> [--perturb INDEX] [--no-write]
"""
import json
import math
import os
import struct
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import K, MOVES, MIDX, DRIFT, SCALE, heights, windows, legal_moves_window, code, decode, apply_move
import law as L

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("JB_WORK", os.path.join(HERE, "work"))


def read_parent(path):
    f = open(path, "rb")
    assert f.read(8) == b"IWPAR001"
    Kf, r, count, dlen = struct.unpack("<IIII", f.read(16))
    D = int.from_bytes(f.read(dlen), "little")
    recs = []
    for _ in range(count):
        level, move, sign, zero, cd, nb = struct.unpack("<BBBBII", f.read(12))
        assert zero == 0 and sign in (0, 1) and move < 6
        N = int.from_bytes(f.read(nb), "little")
        recs.append((level, move, -N if sign else N, cd))
    assert f.read(1) == b"", "trailing bytes"
    return Kf, r, D, recs


class Ratios:
    """beta(hs) / beta(hs[a : len-b]) exactly, from the memory-three cylinder formula, cached on the relevant heights."""
    def __init__(self, law):
        self.q, self.p = law["q"], law["p"]
        self.cache = {}

    def t(self, hs5):
        return self.p[hs5] / self.q[hs5[:4]]

    def ratio(self, hs, a, b):
        k = len(hs) - 1
        key = (tuple(hs[:a + 4]), tuple(hs[k - b - 3:]), a, b, k)
        v = self.cache.get(key)
        if v is None:
            v = self.q[tuple(hs[0:4])] / self.q[tuple(hs[a:a + 4])]
            for i in range(4, a + 4):
                v *= self.t(tuple(hs[i - 4:i + 1]))
            for i in range(k - b + 1, k + 1):
                v *= self.t(tuple(hs[i - 4:i + 1]))
            self.cache[key] = v
        return v


def main():
    args = sys.argv[1:]
    perturb = None
    if "--perturb" in args:
        perturb = int(args[args.index("--perturb") + 1])
    t0 = time.time()
    Kf, r, D, recs = read_parent(args[0])
    assert Kf == K and r == 6
    n = 2 * r + 1
    print("records", len(recs), "D bits", D.bit_length())
    if perturb is not None:
        lv, mv, N, cd = recs[perturb]
        recs[perturb] = (lv, mv, N + 1, cd)
        print("PERTURBED record", perturb, "numerator +1")
    # C3: canonical order and exactly the legal pairs
    expected = []
    for (lev, w) in windows(n):
        for mv in legal_moves_window(lev, w):
            expected.append((lev, MIDX[mv], code(w)))
    got = [(lv, mv, cd) for (lv, mv, N, cd) in recs]
    order_ok = got == expected
    print("windows", len(windows(n)), "legal pairs", len(expected), "records equal legal pairs in canonical order:", order_ok)
    nonzero = all(N != 0 for (_, _, N, _) in recs)
    neg = sum(1 for (_, _, N, _) in recs if N < 0)
    sumN = sum(N for (_, _, N, _) in recs)
    da = sum(N * DRIFT[MOVES[mv]][0] for (_, mv, N, _) in recs)
    db = sum(N * DRIFT[MOVES[mv]][1] for (_, mv, N, _) in recs)
    print("C3: sum N == D:", sumN == D, " drift a == 0:", da == 0, " drift b == 0:", db == 0,
          " negative:", neg, " none zero:", nonzero)
    J = {}
    for (lv, mv, N, cd) in recs:
        J[(lv, cd, mv)] = N
    law = L.build_law()
    R = Ratios(law)
    # C2: rows (S_6). Multiply SPEC 7.1 by D*beta(Y):
    #   N(Y_L,L) b(Y)/b(Y_L) + N(Y_R,R) b(Y)/b(Y_R) + N(w_S,S) + sum N(w_M,M) b(V_M)/b(w_M) - sum_xi N(u,xi) b(Y)/b(u) = 0
    rows = 0; bad = 0; badlist = []
    for (lev, Y) in windows(n + 2):
        rows += 1
        hY = heights(lev, Y)
        terms = []
        # L: window Y_L = Y[2:], level lev + y0 + y1, move L
        lvL = hY[2]
        NL = J.get((lvL, code(Y[2:]), 0))
        if NL is not None:
            terms.append((NL, R.ratio(hY, 2, 0)))
        NR = J.get((lev, code(Y[:n]), 1))
        if NR is not None:
            terms.append((NR, R.ratio(hY, 0, 2)))
        mid = Y[r:r + 3]
        for pat, rep, mv in (("+-+", "+", "S1"), ("+--", "-", "S2"), ("--+", "-", "S3")):
            if mid == pat:
                wS = Y[:r] + rep + Y[r + 3:]
                if apply_move(lev, wS, r, mv) is not None:
                    NS = J.get((lev, code(wS), MIDX[mv]))
                    assert NS is not None
                    terms.append((NS, Fr(1)))
        y = Y[r + 1]
        pairs = [("-", "-")] if y == "-" else [("-", "+"), ("+", "-")]
        for (pp, qq) in pairs:
            V = Y[:r + 1] + pp + "+" + qq + Y[r + 2:]
            hV = heights(lev, V)
            if all(0 <= h < K for h in hV):
                wM = V[2:n + 2]
                NM = J.get((hV[2], code(wM), 5))
                assert NM is not None, "merge source missing"
                terms.append((NM, R.ratio(hV, 2, 2)))
        u = Y[1:n + 1]
        ru = R.ratio(hY, 1, 1)
        for mv in range(6):
            Nu = J.get((hY[1], code(u), mv))
            if Nu is not None:
                terms.append((-Nu, ru))
        Lc = 1
        for (_, f) in terms:
            Lc = Lc * f.denominator // math.gcd(Lc, f.denominator)
        tot = sum(Nv * (f.numerator * (Lc // f.denominator)) for (Nv, f) in terms)
        if tot != 0:
            bad += 1
            if len(badlist) < 5:
                badlist.append((lev, Y))
    print(f"C2: rows {rows}, rows with nonzero residual {bad}", badlist)
    print("time", round(time.time() - t0, 1), "s")
    rep = dict(records=len(recs), D_bits=D.bit_length(), order_ok=order_ok, sumN_eq_D=sumN == D, drift_a=da == 0,
               drift_b=db == 0, negative=neg, nonzero=nonzero, rows=rows, bad_rows=bad, perturbed=perturb)
    if "--no-write" not in args:
        primes = json.load(open(os.path.join(WORK, "primes.json")))["primes"]
        outname = args[args.index("--out") + 1] if "--out" in args else "parent_c.bin"
        out = open(os.path.join(WORK, outname), "wb")
        out.write(struct.pack("<I", len(recs)))
        bcache = {}
        for (lv, mv, N, cd) in recs:
            w = decode(cd, n)
            b = bcache.get((lv, cd))
            if b is None:
                b = L.beta_heights(law, heights(lv, w)); bcache[(lv, cd)] = b
            # x_P = N / (D * b)
            num = N * b.denominator
            den = D * b.numerator
            lo = (num << SCALE) // den
            mods = []
            for P in primes:
                mods.append((num % P) * pow(den % P, -1, P) % P * (1 << 64) % P)
            out.write(struct.pack("<BBHI", lv, mv, 0, cd))
            out.write((lo & ((1 << 128) - 1)).to_bytes(16, "little"))
            out.write(struct.pack("<QQ", mods[0], mods[1]))
        out.close()
        print("wrote work/" + outname)
    ok = order_ok and sumN == D and da == 0 and db == 0 and bad == 0 and nonzero
    print("PARENT", "PASS" if ok else "FAIL", json.dumps(rep))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
