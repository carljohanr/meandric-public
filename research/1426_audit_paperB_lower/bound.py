"""1426 (JB): turn a flowcheck log into a bound, with exact rational comparisons and own outward-rounded logs.

h_lo = H_lo / M_hi (both exact dyadic rationals printed by flowcheck in Q120). The claim e^{2 h_lo} >= t is checked as
2*h_lo >= U(t), U(t) an exact rational upper bound of ln t (atanh series with a rational tail bound).
Also compares the mass interval and the mass mod p with the closed form M = 1 + sum c_d k_d beta(A_d) (walks.py).

Usage: python3 bound.py <flowcheck log> <stage index> [target as fraction or decimal] [walks summary json]
"""
import json
import sys
from fractions import Fraction as Fr

Q = 1 << 120


def atanh_bounds(z, J=60):
    """exact rational lo <= atanh(z) <= hi for 0 <= z < 1/2"""
    s = Fr(0); pw = z
    for j in range(J):
        s += pw / (2 * j + 1); pw *= z * z
    tail = pw / (2 * J + 1) / (1 - z * z)
    return s, s + tail


def ln_bounds(x):
    """exact rational bounds of ln x for x > 0"""
    x = Fr(x)
    l2lo, l2hi = [2 * v for v in atanh_bounds(Fr(1, 3))]
    k = 0
    while x > Fr(4, 3): x /= 2; k += 1
    while x < Fr(2, 3): x *= 2; k -= 1
    z = (x - 1) / (x + 1)
    if z >= 0:
        lo, hi = atanh_bounds(z)
    else:
        lo, hi = atanh_bounds(-z); lo, hi = -hi, -lo
    lo, hi = 2 * lo, 2 * hi
    if k >= 0: return lo + k * l2lo, hi + k * l2hi
    return lo + k * l2hi, hi + k * l2lo


def exp_lo(y, terms=80):
    """exact rational lower bound of e^y for y > 0 (partial Taylor sum)"""
    s = Fr(0); t = Fr(1)
    for i in range(terms):
        s += t; t = t * y / (i + 1)
    return s


def dec_floor(x, digits):
    """decimal string of x >= 0 truncated (rounded down) to the given number of digits after the point"""
    v = x.numerator * 10**digits // x.denominator
    s = str(v).rjust(digits + 1, "0")
    return s[:-digits] + "." + s[-digits:]


def parse(logf):
    d = {}
    for line in open(logf):
        p = line.split()
        if not p: continue
        if p[0] in ("H_lo_q120", "M_hi_q120", "M_lo_q120"): d[p[0]] = int(p[1])
        if p[0] == "mass_mod_montgomery": d["mass_mod"] = [int(p[1]), int(p[2])]
        if p[0] == "primes": d["primes"] = [int(p[1]), int(p[2])]
        if line.startswith("pass1"): d["pass1"] = line.strip()
        if line.startswith("pass2"): d["pass2"] = line.strip()
        if line.startswith("drift_mod"): d["drift"] = line.strip()
    return d


def main():
    logf, stage = sys.argv[1], int(sys.argv[2])
    target = Fr(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "-" else None
    wsum = sys.argv[4] if len(sys.argv) > 4 else "work/edges_c_summary.json"
    d = parse(logf)
    H, Mhi, Mlo = d["H_lo_q120"], d["M_hi_q120"], d["M_lo_q120"]
    h = Fr(H, Mhi)
    print("h_lo = H_lo/M_hi =", dec_floor(h, 40))
    print("M in [", float(Fr(Mlo, Q)), ",", float(Fr(Mhi, Q)), "], width", float(Fr(Mhi - Mlo, Q)))
    W = json.load(open(wsum))
    if stage >= 0:
        ps = W["per_stage"][stage]
        cf_lo = Fr(int(ps["M_lo"]), 1 << W["scale"]); cf_hi = Fr(int(ps["M_hi"]), 1 << W["scale"])
        Mmod = ps["M_mod"]
    else:
        cf_lo = cf_hi = Fr(1); Mmod = [1, 1]
    print("closed-form mass inside [M_lo, M_hi]:", Fr(Mlo, Q) <= cf_lo and cf_hi <= Fr(Mhi, Q))
    mm = [(m * pow(1 << 64, -1, P)) % P for m, P in zip(d["mass_mod"], d["primes"])]
    print("mass mod p equals closed form mod p:", mm == Mmod)
    two_h = 2 * h
    e2h = exp_lo(two_h)
    print("e^{2 h_lo} >= (Taylor lower bound)", dec_floor(e2h, 15))
    if target is not None:
        lo, hi = ln_bounds(target)
        ok = two_h >= hi
        slack = float(h - hi / 2)
        print(f"target {target} = {float(target)!r}: 2 h_lo >= upper bound of ln(target): {ok}; slack in h >= {slack:.4e}")
        return ok
    return True


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
