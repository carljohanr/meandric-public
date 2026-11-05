"""1390 (HP): second, literal check of Paper E Lemma 5 on random states, in Python Fractions.

v(x) is computed from scratch on the whole stacks (not as a ratio), and sum_y v(y) <= R v(x) is tested
exactly.  Also: z_1..z_10 by the transfer of Paper E Section 2, z_n <= R^(2n) for those n, and R^2 exactly.
Written from Paper E and the BICCERT1 format only.
Usage: python hp_upperE_literal.py cert.bin [nsamples]
"""
import sys, struct, random
from fractions import Fraction as F
import numpy as np

path = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
with open(path, 'rb') as f:
    head = f.read(84)
assert head[:8] == b'BICCERT1'
r = struct.unpack('<i', head[8:12])[0]
an, ad, pn, pd, tn, td, P, Q, cnt = struct.unpack('<9Q', head[12:84])
tab = np.memmap(path, dtype='<u4', mode='r', offset=84, shape=(cnt,))
a, p, t, R = F(an, ad), F(pn, pd), F(tn, td), F(P, Q)

def view(S):  # S: tuple bottom-first
    if len(S) < r:
        return ('^',) + S
    return S[len(S) - r:]

def code(V):
    s = '1' + ''.join(str(x) for x in V if x != '^')
    return int(s, 2)  # leading 1, then parities bottom-to-top, so the top is the lowest bit

def u(VU, VL, c):
    if c == 1:
        VU = tuple(x if x == '^' else 1 - x for x in VU)
        VL = tuple(x if x == '^' else 1 - x for x in VL)
    return int(tab[(code(VU) << (r + 1)) | code(VL)])

def E(S):
    return sum(1 for i in range(len(S) - 1) if S[i] == S[i + 1])

def T(S):
    return sum(1 for i in range(len(S) - 2) if S[i] == S[i + 1] == S[i + 2])

def v(U, L, c):
    return a ** (len(U) + len(L)) * p ** (E(U) + E(L)) * t ** (T(U) + T(L)) * u(view(U), view(L), c)

def succ(U, L, c):
    out = [(U + (c,), L, 1 - c), (U, L + (c,), 1 - c)]
    if U and U[-1] == 1 - c:
        out.append((U[:-1], L, 1 - c))
    if L and L[-1] == 1 - c:
        out.append((U, L[:-1], 1 - c))
    return out

# z_n by the transfer
from collections import Counter
cur = Counter({((), (), 0): 1})
z = []
for step in range(1, 21):
    nxt = Counter()
    for (U, L, c), m in cur.items():
        for y in succ(U, L, c):
            nxt[y] += m
    cur = nxt
    if step % 2 == 0:
        z.append(cur[((), (), 0)])
print('z_1..z_10 =', z)
assert z == [2, 8, 40, 228, 1424, 9520, 67064, 492292, 3735112, 29114128]
print('z_n <= R^(2n) for n<=10:', all(F(zn) <= R ** (2 * n) for n, zn in enumerate(z, 1)))
R2 = R * R
print('R^2 =', R2, '=', f"{R2.numerator * 10**21 // R2.denominator}e-21 (floor)")

rnd = random.Random(1390)
worst = F(0); bad = 0
for i in range(N):
    mode = i % 3
    if mode == 0:   # arbitrary state, long stacks
        U = tuple(rnd.randint(0, 1) for _ in range(rnd.randint(0, 30)))
        L = tuple(rnd.randint(0, 1) for _ in range(rnd.randint(0, 30)))
        c = rnd.randint(0, 1)
    elif mode == 1:  # charge-consistent state reached by a random walk from e
        U, L, c = (), (), 0
        for _ in range(rnd.randint(1, 80)):
            U, L, c = rnd.choice(succ(U, L, c))
    else:            # near-tight family: alternating stacks with runs of equal symbols
        U = tuple(rnd.choice([0, 1]) for _ in range(rnd.randint(0, 3)))
        L = tuple([1] * rnd.randint(0, 6) + [0] * rnd.randint(0, 8))
        c = 0
    vx = v(U, L, c)
    q = sum(v(*y) for y in succ(U, L, c)) / vx
    if q > R:
        bad += 1
    if q > worst:
        worst = q
print(f'random literal states: {N}, violations: {bad}, largest quotient seen: {float(worst):.12f}')
