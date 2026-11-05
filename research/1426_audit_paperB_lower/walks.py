"""1426 (JB): replay every closed-walk direction (check C4) and write the edge terms for the C checker.

For each direction d of stage s (radius r_s): the walk is replayed from its anchor; asserted: every move legal
(interior head, move applicable, new word in the strip), closed, protected (first 3, last 3 letters and letter sum
unchanged at every tape), every source tape fits r_s <= 10, zero drift. Coefficient columns must agree exactly.
Edge term: c_d * beta(A_d)/beta(T_{j-1}) (exact rational), written as floor(term * 2^88) and term mod p (Montgomery).
Also accumulates the closed-form mass M = 1 + sum_d c_d k_d beta(A_d) as an exact fixed-point interval
(scale 2^-256) and mod p.

Usage: python3 walks.py <export dir> [--flip STAGE INDEX] [--no-write]
Output: work/edges_c.bin (header u64 count; records of 44 bytes), logs on stdout.
"""
import json
import os
import struct
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import K, MOVES, MIDX, DRIFT, SCALE, heights, in_strip, apply_move, code
import law as L

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("JB_WORK", os.path.join(HERE, "work"))
FILES = ["dirs_0_base_r8.tsv", "dirs_1_sparse_r9.tsv", "dirs_2_aggregated_r10.tsv", "dirs_3_increment1_r10.tsv",
         "dirs_4_increment2_r10.tsv"]
RADII = [8, 9, 10, 10, 10]
MSCALE = 256


def main():
    args = sys.argv[1:]
    exp = args[0]
    flip = None
    if "--flip" in args:
        i = args.index("--flip"); flip = (int(args[i + 1]), int(args[i + 2]))
    write = "--no-write" not in args
    scale = None
    if "--scale-coef" in args:  # damage test: multiply one coefficient by an integer factor
        i = args.index("--scale-coef"); scale = (int(args[i + 1]), int(args[i + 2]), int(args[i + 3]))
    force = "--force-include" in args  # damage test: write the edges of a rejected walk anyway
    outname = args[args.index("--out") + 1] if "--out" in args else "edges_c.bin"
    law = L.build_law()
    q, p = law["q"], law["p"]
    tN = {hs: (p[hs] / q[hs[:4]]).numerator for hs in p}
    tD = {hs: (p[hs] / q[hs[:4]]).denominator for hs in p}
    primes = json.load(open(os.path.join(WORK, "primes.json")))["primes"]
    stages = json.load(open(os.path.join(exp, "stages.json")))
    out = open(os.path.join(WORK, outname), "wb") if write else None
    if out: out.write(struct.pack("<Q", 0))
    nedges = 0
    Mlo = 1 << MSCALE; Mhi = 1 << MSCALE  # mass 1 from the parent (C3)
    Mmod = [1, 1]
    total_dirs = 0; failures = []
    ratio_cache = {}
    per_stage = []
    max_abs_term = 0.0
    t0 = time.time()
    for s, fn in enumerate(FILES):
        rs = RADII[s]
        shapes = {}
        nd = 0
        for line in open(os.path.join(exp, fn)):
            if line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            assert len(cols) == 10
            idx, lev, A, h, mv, chex, cnum, cpow = int(cols[0]), int(cols[1]), cols[2], int(cols[3]), cols[4].split(","), cols[5], int(cols[6]), int(cols[7])
            if flip is not None and flip == (s, idx):
                # damage: replace the first head move by its opposite
                cand = [j for j, m in enumerate(mv) if m in ("L", "R")]
                j = cand[-1] if "--flip-last" in args else cand[0]
                mv = list(mv); mv[j] = "R" if mv[j] == "L" else "L"
                print("FLIPPED", s, idx, "move", j, "->", ",".join(mv))
            if scale is not None and (s, idx) == scale[:2]:
                e2 = scale[2]  # multiply by 2^e2 (exact in float64)
                if cpow >= e2: cpow -= e2
                else: cnum <<= (e2 - cpow); cpow = 0
                chex = float(Fr(cnum, 1 << cpow)).hex()
                print("SCALED coefficient", s, idx, "by 2 **", e2)
            ok = True; why = ""
            if idx != nd: ok = False; why = "index"
            c = Fr(cnum, 1 << cpow)
            if c == 0 or Fr(float.fromhex(chex)) != c or cpow < 0: ok = False; why = "coefficient"
            # replay
            tapes = [(A, h)]
            if not in_strip(lev, A): ok = False; why = "anchor not in strip"
            w, hh = A, h
            da = db = 0
            for m in mv:
                if m not in MIDX: ok = False; why = "bad move name"; break
                if not (hh <= rs and len(w) - 1 - hh <= rs): ok = False; why = "source does not fit"; break
                nxt = apply_move(lev, w, hh, m)
                if nxt is None: ok = False; why = "illegal move"; break
                w, hh = nxt
                da += DRIFT[m][0]; db += DRIFT[m][1]
                tapes.append((w, hh))
            if ok:
                if (w, hh) != (A, h): ok = False; why = "not closed"
                elif (da, db) != (0, 0): ok = False; why = "drift"
                elif rs > 10: ok = False; why = "radius"
                else:
                    sA = sum(1 if ch == "+" else -1 for ch in A)
                    for (T, _) in tapes:
                        if T[:3] != A[:3] or T[-3:] != A[-3:] or sum(1 if ch == "+" else -1 for ch in T) != sA:
                            ok = False; why = "not protected"; break
            if not ok:
                failures.append((s, idx, why))
                if not (force and why in ("not closed", "drift")):
                    nd += 1
                    continue
                print("FORCE-INCLUDING rejected walk", s, idx, why, "tapes", tapes)
            key = (len(A), h, ",".join(mv)); shapes[key] = shapes.get(key, 0) + 1
            if ok: nd += 1
            else: nd += 1
            hA = heights(lev, A)
            fA = [tuple(hA[i - 4:i + 1]) for i in range(4, len(hA))]
            # beta(A) as NUM/DEN (unreduced), for the mass
            bN = q[tuple(hA[:4])].numerator; bD = q[tuple(hA[:4])].denominator
            for f in fA:
                bN *= tN[f]; bD *= tD[f]
            k = len(mv)
            num = cnum * k * bN; den = bD << cpow
            Mlo += (num << MSCALE) // den
            Mhi += -((-(num << MSCALE)) // den)
            for i, P in enumerate(primes):
                Mmod[i] = (Mmod[i] + (num % P) * pow(den % P, -1, P)) % P
            for j, m in enumerate(mv[:len(tapes) - 1] if not ok else mv):
                T, hT = tapes[j]
                if T == A:
                    NUM = DEN = 1
                else:
                    hT_ = heights(lev, T)
                    fT = [tuple(hT_[i - 4:i + 1]) for i in range(4, len(hT_))]
                    a0 = 0
                    while a0 < min(len(fA), len(fT)) and fA[a0] == fT[a0]: a0 += 1
                    b0 = 0
                    while b0 < min(len(fA), len(fT)) - a0 and fA[-1 - b0] == fT[-1 - b0]: b0 += 1
                    ka = tuple(fA[a0:len(fA) - b0]); kt = tuple(fT[a0:len(fT) - b0])
                    if hA[:4] != hT_[:4]: raise AssertionError("first quartet differs")
                    rc = ratio_cache.get((ka, kt))
                    if rc is None:
                        NUM = 1; DEN = 1
                        for f in ka: NUM *= tN[f]; DEN *= tD[f]
                        for f in kt: NUM *= tD[f]; DEN *= tN[f]
                        rc = (NUM, DEN); ratio_cache[(ka, kt)] = rc
                    NUM, DEN = rc
                tn = cnum * NUM; td = DEN << cpow
                lo = (tn << SCALE) // td
                max_abs_term = max(max_abs_term, abs(lo) / 2.0 ** SCALE)
                mods = [(tn % P) * pow(td % P, -1, P) % P * (1 << 64) % P for P in primes]
                if out:
                    out.write(struct.pack("<BBBBBBHI", rs, lev, len(T), hT, MIDX[m], s, 0, code(T)))
                    out.write((lo & ((1 << 128) - 1)).to_bytes(16, "little"))
                    out.write(struct.pack("<QQ", mods[0], mods[1]))
                nedges += 1
        st = stages["stages"][s]
        exp_shapes = {(wt["anchor_length"], wt["head"], wt["moves"]): wt["count"] for wt in st["walk_types"]}
        print(f"stage {s} {fn}: directions {nd} (stages.json {st['directions']}), shapes {len(shapes)}, "
              f"shapes equal stages.json: {shapes == exp_shapes}, time {time.time()-t0:.0f}s", flush=True)
        total_dirs += nd
        per_stage.append({"stage": s, "M_lo": str(Mlo), "M_hi": str(Mhi), "M_mod": list(Mmod)})
    if out:
        out.seek(0); out.write(struct.pack("<Q", nedges)); out.close()
    print("directions", total_dirs, "edges", nedges, "failures", len(failures), failures[:10])
    print("max |edge density term|", max_abs_term)
    print("closed-form mass M in [", Mlo / 2.0 ** MSCALE, "]; exact interval numerators at scale 2^-256:")
    print("  M_lo", Mlo)
    print("  M_hi", Mhi)
    from decimal import Decimal, getcontext
    getcontext().prec = 80
    print("  M ~", Decimal(Mlo) / Decimal(2) ** MSCALE)
    print("  M mod p (plain form)", Mmod)
    if write: json.dump({"M_lo": str(Mlo), "M_hi": str(Mhi), "scale": MSCALE, "M_mod": Mmod, "edges": nedges,
               "directions": total_dirs, "failures": failures[:100], "per_stage": per_stage}, open(os.path.join(WORK, outname.replace(".bin", "_summary.json")), "w"))
    sys.exit(0 if not failures else 1)


if __name__ == "__main__":
    main()
