"""1426 (JB): test the C primitives against exact values (run with the venv python, which has mpmath).

- log_ratio_lo(A, B): must satisfy L/2^64 <= ln(A/B) and be within 2^-56 of it, on random and extreme inputs;
- beta_mant: lo*2^e <= beta(word) <= hi*2^e with relative width < 1e-16, on random words of 13 to 25 letters;
- beta_mod: equals beta(word) mod p (Montgomery form).
Usage: python primtest.py
"""
import json
import os
import random
import subprocess
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, log
import law as L
from common import heights, in_strip

mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(1426)
lines = []; logcases = []; betacases = []
for i in range(20000):
    kind = i % 4
    if kind == 0:
        B = rng.randrange(1, 1 << rng.randrange(1, 127)); A = B + rng.randrange(0, B + 1)
    elif kind == 1:
        B = rng.randrange(1 << 80, 1 << 100); A = B * rng.randrange(1, 10**6) + rng.randrange(0, B)
    elif kind == 2:
        B = rng.randrange(1 << 82, 1 << 92); A = B + rng.randrange(0, 1 << rng.randrange(1, 60))
    else:
        A = rng.randrange(1, 1 << 127); B = rng.randrange(1, 1 << 127)
    if A >= 1 << 128 or B >= 1 << 128: continue
    logcases.append((A, B))
    lines.append(f"L {A >> 64} {A & ((1 << 64) - 1)} {B >> 64} {B & ((1 << 64) - 1)}")
for i in range(5000):
    n = rng.randrange(13, 26)
    while True:
        lev = rng.randrange(12); w = "".join(rng.choice("+-") for _ in range(n))
        if in_strip(lev, w): break
    c = int(w.replace("+", "1").replace("-", "0"), 2)
    betacases.append((lev, w)); lines.append(f"B {lev} {c} {n}")
out = subprocess.run([os.path.join(HERE, "work", "primtest"), os.path.join(HERE, "work", "law_c.txt")],
                     input="\n".join(lines) + "\n", capture_output=True, text=True, check=True).stdout.split("\n")
res_l = [l for l in out if l.startswith("L ")]; res_b = [l for l in out if l.startswith("B ")]
law = L.build_law()
primes = json.load(open(os.path.join(HERE, "work", "primes.json")))["primes"]
bad = 0; worst = mpf(0)
for (A, B), l in zip(logcases, res_l):
    Lq = int(l.split()[1])
    true = log(mpf(A) / mpf(B)) if A > B else mpf(0)
    got = mpf(Lq) / mpf(2) ** 64
    if got > true: bad += 1
    worst = max(worst, true - got)
print("log cases", len(logcases), "violations (bound above true value)", bad, "largest gap", float(worst))
badb = 0; worstw = 0.0; badm = 0
for (lev, w), l in zip(betacases, res_b):
    p = l.split(); lo, elo, hi, ehi, m0, m1 = int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6])
    b = L.beta_heights(law, heights(lev, w))
    blo = Fr(lo) * Fr(2) ** elo; bhi = Fr(hi) * Fr(2) ** ehi
    if not (blo <= b <= bhi): badb += 1
    worstw = max(worstw, float((bhi - blo) / b))
    for m, P in zip((m0, m1), primes):
        if m != (b.numerator * pow(b.denominator, -1, P) % P) * (1 << 64) % P: badm += 1
print("beta cases", len(betacases), "enclosure violations", badb, "largest relative width", worstw, "mod p mismatches", badm)
sys.exit(1 if bad or badb or badm or worst > mpf(2) ** -56 else 0)
