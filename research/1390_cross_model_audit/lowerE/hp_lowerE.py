"""1390 (HP): independent checker for Paper E, Theorem L and Proposition 7 (reversible flows),
written from the paper (Section 4 and the appendix "Format of the lower certificates") only.
No project checker code was read.

Steps
 1. read BICFLOW1, check sizes, positivity of every X;
 2. build the flow J of Proposition 7 in exact integers (common scale 8 sd^(2r-1) dropped);
 3. check (P), (D) exactly, and (S) literally WITHOUT using Proposition 7: for every window, every
    allowed move and every value of the four hidden symbols (two below the window on each stack, drawn
    from the reference chain), the move is applied to the concrete top r+2 symbols, and the probability
    of each extended window (r+1 symbols per stack) after the move is compared with its probability before;
 4. bound h from below with an integer fixed-point atanh series (every rounding downward), and compare
    e^{2h} with the target using exact rationals and rigorous series enclosures.
Usage: python hp_lowerE.py flow.bin target [--noS] [--damage MODE]
"""
import sys, struct, time
from fractions import Fraction as F
import numpy as np

t0 = time.time()
path = sys.argv[1]
target = F(sys.argv[2])
doS = '--noS' not in sys.argv
damage = sys.argv[sys.argv.index('--damage') + 1] if '--damage' in sys.argv else None

raw = open(path, 'rb').read()
def reject(msg):
    print('REJECT:', msg); sys.exit(1)
if raw[:8] != b'BICFLOW1': reject('bad magic')
r = struct.unpack('<i', raw[8:12])[0]
if not 1 <= r <= 12: reject('radius out of range')
sn, sd, nA, nB = struct.unpack('<4q', raw[12:44])
print(f'r={r} s={sn}/{sd} nA={nA} nB={nB}')
if nA != 4 ** r or nB != 4 ** r: reject('block sizes != 4^r')
if not (0 < sn < sd): reject('s not in (0,1)')
if len(raw) != 44 + 8 * (nA + nB): reject(f'file length {len(raw)} != {44 + 8 * (nA + nB)}')
X = np.frombuffer(raw, dtype='<i8', offset=44, count=nA + nB)
if damage == 'xzero':
    X = X.copy(); X[12345] = 0
if (X <= 0).any(): reject(f'{int((X <= 0).sum())} nonpositive X')
XA = X[:nA].astype(object); XB = X[nA:].astype(object)
print('X: all positive, max', int(X.max()), 'bits', int(X.max()).bit_length())

M = (1 << r) - 1
NW = 1 << (2 * r + 1)
idx = np.arange(NW, dtype=np.int64)
U = idx & M; L = (idx >> r) & M; c = idx >> (2 * r)
U1 = U & 1; L1 = L & 1
sk = sd - sn  # numerator of 1-s (probability of no switch)

def switches(w):  # number of adjacent unequal pairs in an r-letter word
    return np.bitwise_count((w ^ (w >> 1)) & ((1 << (r - 1)) - 1)).astype(np.int64)

sw = switches(U) + switches(L)
PW = np.array([sn ** k * sk ** (2 * r - 2 - k) for k in range(2 * r - 1)], dtype=object)
Bw = PW[sw]           # beta(w) * 8 sd^(2r-2)
H = 1 << (2 * r - 1)  # c * 2^(2r-1) offset in the X blocks

# window after pushing c on U: U' = (c, U_1..U_{r-1}), parity 1-c; its X_U index uses U'_2..U'_r = U_1..U_{r-1}
iA_push = (1 - c) * H + L * (1 << (r - 1)) + (U & ((1 << (r - 1)) - 1))
iA_pop = c * H + L * (1 << (r - 1)) + (U >> 1)
iB_push = (1 - c) * H + (L & ((1 << (r - 1)) - 1)) * (1 << r) + U
iB_pop = c * H + (L >> 1) * (1 << r) + U
qU = np.where(c != U1, sn, sk).astype(object)   # q(c | w) * sd
qL = np.where(c != L1, sn, sk).astype(object)
popA = (U1 == 1 - c); popB = (L1 == 1 - c)
zero = np.zeros(NW, dtype=object); zero[:] = 0

J = {}
J['pA'] = Bw * qU * XA[iA_push]
J['pB'] = Bw * qL * XB[iB_push]
J['qA'] = np.where(popA, Bw * sd * XA[iA_pop], zero)
J['qB'] = np.where(popB, Bw * sd * XB[iB_pop], zero)
if damage == 'noq':      # pushes without the reference factor q
    J['pA'] = Bw * sd * XA[iA_push]
if damage == 'one':      # one flow value changed by one unit
    J['pA'][777] = J['pA'][777] + 1
nu = J['pA'] + J['pB'] + J['qA'] + J['qB']
Z = int(nu.sum())
nflow = NW * 2 + int(popA.sum()) + int(popB.sum())
print(f'windows {NW}, flow values {nflow}; built in {time.time() - t0:.1f}s')

# (P) and (D)
assert all(int(v) > 0 for v in J['pA']) and all(int(v) > 0 for v in J['pB']), '(P) fails'
print('(P): all push flows positive')
dA = int(J['pA'].sum()) - int(J['qA'].sum()); dB = int(J['pB'].sum()) - int(J['qB'].sum())
print(f'(D): sum pA - sum qA = {dA}, sum pB - sum qB = {dB}')
okD = (dA == 0 and dB == 0)

# (S), literally
okS = None
if doS:
    t1 = time.time()
    M1 = (1 << (r + 1)) - 1
    NY = 1 << (2 * r + 3)
    before = np.zeros(NY, dtype=object); before[:] = 0
    after = np.zeros(NY, dtype=object); after[:] = 0
    Ur = (U >> (r - 1)) & 1; Lr = (L >> (r - 1)) & 1
    def kap(x, y):  # numerator (times sd) of P(next = x | previous = y)
        return np.where(x != y, sn, sk).astype(object)
    revealed_pop = 'one' if damage == 'onefresh' else 'two'
    for h1 in (0, 1):
        for h2 in (0, 1):
            for l1 in (0, 1):
                for l2 in (0, 1):
                    wt = kap(h1, Ur) * kap(h2, np.full(NW, h1)) * kap(l1, Lr) * kap(l2, np.full(NW, l1))
                    Uf = U | (h1 << r) | (h2 << (r + 1))
                    Lf = L | (l1 << r) | (l2 << (r + 1))
                    def Yi(Ue, Le, cc):
                        return cc * (1 << (2 * r + 2)) + Le * (1 << (r + 1)) + Ue
                    # before: the present extended window
                    Y = Yi(Uf & M1, Lf & M1, c)
                    np.add.at(before, Y, nu * wt)
                    # pushes
                    Y = Yi(((Uf << 1) | c) & M1, Lf & M1, 1 - c)
                    np.add.at(after, Y, J['pA'] * wt)
                    Y = Yi(Uf & M1, ((Lf << 1) | c) & M1, 1 - c)
                    np.add.at(after, Y, J['pB'] * wt)
                    # pops (only where allowed)
                    if revealed_pop == 'two':
                        Ue = (Uf >> 1) & M1
                    else:  # wrong model: one fresh symbol, the second hidden symbol forced to equal the first
                        Ue = ((Uf >> 1) & M1)
                        Ue = np.where(h2 == h1, Ue, -1)
                    sel = popA if revealed_pop == 'two' else (popA & (Ue >= 0))
                    np.add.at(after, Yi(Ue[sel], Lf[sel] & M1, (1 - c)[sel]), (J['qA'] * wt)[sel])
                    Le = (Lf >> 1) & M1
                    np.add.at(after, Yi(Uf[popB] & M1, Le[popB], (1 - c)[popB]), (J['qB'] * wt)[popB])
    diff = before - after
    nbad = int(sum(1 for v in diff if v != 0))
    okS = (nbad == 0)
    print(f'(S): {NY} extended-window equations checked exactly, {nbad} fail; {time.time() - t1:.1f}s')

# entropy lower bound
P = 96
ONE = 1 << P
def ln2_bounds():
    # ln 2 = 2 atanh(1/3); partial sums are lower bounds, tail <= z^(2K+3)/((2K+3)(1-z^2))
    z = F(1, 3); s = F(0); K = 80
    for j in range(K + 1):
        s += z ** (2 * j + 1) / (2 * j + 1)
    tail = z ** (2 * K + 3) / ((2 * K + 3) * (1 - z * z))
    return 2 * s, 2 * (s + tail)
LN2_lo, LN2_hi = ln2_bounds()
LN2_lo_fx = (LN2_lo.numerator * ONE) // LN2_lo.denominator

def log_lo_fx(n, d):
    """floor-type fixed-point lower bound (scaled by 2^P) for log(n/d), n >= d > 0."""
    k = n.bit_length() - d.bit_length()
    if k < 0 or (d << k) > n:
        k -= 1
    if k < 0:
        raise ValueError('ratio < 1')
    dk = d << k
    z = ((n - dk) << P) // (n + dk)
    z2 = (z * z) >> P
    acc = z; term = z; j = 1
    while term:
        term = (term * z2) >> P
        acc += term // (2 * j + 1)
        j += 1
    return k * LN2_lo_fx + 2 * acc

t2 = time.time()
S = 0
for key in ('pA', 'pB', 'qA', 'qB'):
    arr = J[key]
    for jv, nv in zip(arr, nu):
        if jv:
            S += jv * log_lo_fx(nv, jv)
h_lo = F(S, Z * ONE)          # h >= h_lo
print(f'entropy: h >= {float(h_lo):.15f} ({time.time() - t2:.1f}s)')

def exp_lo(x, K=80):  # x > 0: partial sums of the exponential series are lower bounds
    s = F(0); term = F(1)
    for k in range(K):
        s += term; term = term * x / (k + 1)
    return s
def log_hi(y):  # upper bound for log y, y >= 1 rational
    k = 0
    while y >= 2:
        y /= 2; k += 1
    z = (y - 1) / (y + 1); s = F(0); K = 60
    for j in range(K + 1):
        s += z ** (2 * j + 1) / (2 * j + 1)
    tail = z ** (2 * K + 3) / ((2 * K + 3) * (1 - z * z))
    return k * LN2_hi + 2 * (s + tail)

# compare: e^{2h} >= target  <=  2 h_lo >= log_hi(target)
lt = log_hi(target)
ok_target = (2 * h_lo >= lt)
# report a decimal lower bound for e^{2h} (12 decimals, floored); use a coarse rational for speed
h_lo_c = F((h_lo.numerator * 10 ** 30) // h_lo.denominator, 10 ** 30)
e2 = exp_lo(2 * h_lo_c)
print(f'e^(2h) >= {e2.numerator * 10**13 // e2.denominator / 1e13:.13f} (floored, 13 decimals: {e2.numerator * 10**13 // e2.denominator})')
print(f'slack in h vs target: {float(h_lo - lt / 2):.3e}')
print(f'target {sys.argv[2]}: {"MET" if ok_target else "NOT MET"}')
allok = okD and (okS is not False) and ok_target
print('RESULT:', 'ACCEPT' if allok else 'REJECT', f'(S {"checked" if doS else "not checked"}) total {time.time() - t0:.1f}s')
sys.exit(0 if allok else 1)
