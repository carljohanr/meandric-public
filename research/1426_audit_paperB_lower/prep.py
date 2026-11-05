"""1426 (JB): choose two random 62-bit primes and write the law tables for the C checker.

work/primes.json : the primes (chosen with the OS random generator after the export hashes were verified).
work/law_c.txt   : for each quartet (index a*8+steps) and each quintet factor t = p/q (index a*16+steps):
                   value mod p1, mod p2 in Montgomery form (x*2^64 mod p), and 64-bit mantissas lo, hi with exponent e,
                   lo*2^e <= value <= hi*2^e, 2^63 <= lo <= hi < 2^64.
Usage: python3 prep.py [--law-variant NAME]
"""
import json
import os
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import K, random_prime
import law as L

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("JB_WORK", os.path.join(HERE, "work"))


def mant(v):
    """Return (lo, hi, e) with lo*2^e <= v <= hi*2^e, lo, hi in [2^63, 2^64) (hi may equal 2^64 -> renormalized)."""
    assert v > 0
    n, d = v.numerator, v.denominator
    e = n.bit_length() - d.bit_length() - 64
    while True:
        # m = v / 2^e
        if e >= 0:
            num, den = n, d << e
        else:
            num, den = n << (-e), d
        lo = num // den
        if lo >= 1 << 64:
            e += 1; continue
        if lo < 1 << 63:
            e -= 1; continue
        break
    hi = lo if lo * den == num else lo + 1
    if hi == 1 << 64:
        hi = 1 << 63; ehi = e + 1
        # represent hi with the same exponent is impossible; store as lo exponent e and hi exponent e+1
        return lo, e, hi, ehi
    return lo, e, hi, e


def steps_code(hs):
    c = 0
    for i in range(1, len(hs)):
        c = (c << 1) | (1 if hs[i] > hs[i - 1] else 0)
    return c


def main():
    os.makedirs(WORK, exist_ok=True)
    pf = os.path.join(WORK, "primes.json")
    if os.path.exists(pf):
        primes = json.load(open(pf))["primes"]
    else:
        primes = [random_prime(62), random_prime(62)]
        json.dump({"primes": primes, "note": "random 62-bit primes from secrets.randbits, Miller-Rabin with 13 bases"},
                  open(pf, "w"), indent=1)
    law = L.build_law()
    q, p = law["q"], law["p"]
    for P in primes:
        for v in list(q.values()) + list(p.values()):
            assert v.numerator % P and v.denominator % P, "prime divides a law entry"
    lines = [" ".join(str(P) for P in primes)]
    for hs, v in sorted(q.items()):
        lo, elo, hi, ehi = mant(v)
        mods = [(v.numerator * pow(v.denominator, -1, P) % P) * (1 << 64) % P for P in primes]
        lines.append(f"Q {hs[0]*8 + steps_code(hs)} {mods[0]} {mods[1]} {lo} {elo} {hi} {ehi}")
    for hs, v in sorted(p.items()):
        t = v / q[hs[:4]]
        lo, elo, hi, ehi = mant(t)
        mods = [(t.numerator * pow(t.denominator, -1, P) % P) * (1 << 64) % P for P in primes]
        lines.append(f"T {hs[0]*16 + steps_code(hs)} {mods[0]} {mods[1]} {lo} {elo} {hi} {ehi}")
    open(os.path.join(WORK, sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else "law_c.txt"), "w").write("\n".join(lines) + "\n")
    print("primes", primes, "quartets", len(q), "quintet factors", len(p))


if __name__ == "__main__":
    main()
