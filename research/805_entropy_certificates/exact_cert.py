#!/usr/bin/env python3
"""805: exact stationary-flow certificate (generator).

Far tape: flip-symmetric order-2 chain, P(+ | - -) = s, P(+ | - +) = t, P(+ | + -) = 1 - t, P(+ | + +) = 1 - s.
1. float optimum of the flow program (flow2.py), made mirror symmetric;
2. exact integer constraint matrix on mirror orbits; a basis by two pivoted QR factorizations;
3. non-basic flows rounded to multiples of 2^-FB; basic flows from the square system A y = b, solved EXACTLY by
   iterative refinement with exact integer residuals followed by rational reconstruction; the result is checked
   exactly (A N = b D) before it is written;
4. output: pickle with r, s, t and the integer numerators of all flows over one common denominator.
The certificate is checked by verify_cert.py, which shares no code with this file or with flow2.py.

usage: exact_cert.py r s_num s_den t_num t_den out.pkl
"""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import sys, math, time, pickle
from fractions import Fraction as Fr
import numpy as np
import scipy.linalg as sl
from flow2 import Flow

FB = 56                                   # bits of the rounded non-basic flows
t00 = time.time()
def log(*a): print(f'[{time.time()-t00:6.0f}s]', *a, flush=True)

def b3_exact(s, t):
    a = {(0, 0): s, (0, 1): t, (1, 0): 1 - t, (1, 1): 1 - s}
    u = t / (2 * (s + t)); v = s / (2 * (s + t)); b2 = {(0, 0): u, (1, 1): u, (0, 1): v, (1, 0): v}
    return [[[b2[x, y] * (a[x, y] if z else 1 - a[x, y]) for z in (0, 1)] for y in (0, 1)] for x in (0, 1)]

def main():
    r = int(sys.argv[1]); s = Fr(int(sys.argv[2]), int(sys.argv[3])); t = Fr(int(sys.argv[4]), int(sys.argv[5])); out = sys.argv[6]
    b3 = b3_exact(s, t)
    fl = Flow(r, b3=[[[float(z) for z in y] for y in x] for x in b3])
    h, J, nu, rn, it = fl.solve(verbose=False)
    log(f'float optimum: exp(2h) = {math.exp(2*h):.9f}, residual {rn:.1e}, min J {J.min():.2e}')
    fx = Flow(r, b3=b3)                                        # exact transition entries
    n = fx.n; nv = fx.nv
    def rev(code, length): return int(format(code, f'0{length}b')[::-1], 2)
    swap = {'L': 'R', 'R': 'L', 'S1': 'S1', 'S2': 'S3', 'S3': 'S2', 'M': 'M'}
    vid = {(int(w), nm): v for v, (w, nm) in enumerate(zip(fx.vw, fx.vname))}
    mir = [vid[(rev(int(w), n), swap[nm])] for w, nm in zip(fx.vw, fx.vname)]
    orb = [-1] * nv; no = 0
    for v in range(nv):
        if orb[v] < 0: orb[v] = orb[mir[v]] = no; no += 1
    orb = np.array(orb); size = np.bincount(orb, minlength=no)
    J = 0.5 * (J + J[mir]); xstar = np.zeros(no); xstar[orb] = J
    # exact rows on orbits, scaled to integers
    SC = 1
    LT, RT = fx.LT, fx.RT
    for f in (0, 1):
        for x in (0, 1):
            for y in (0, 1): SC = math.lcm(SC, Fr(LT[f][x][y]).denominator, Fr(RT[x][y][f]).denominator)
    SC = SC ** 4
    rows = {}
    ri, ci, val = fx.tri
    keep = [W for W in range(fx.nrow) if W <= rev(W, n + 2)]; kid = {W: k for k, W in enumerate(keep)}
    R = [dict() for _ in keep]
    for row, v, pr in zip(ri, ci, val):
        if row in kid:
            d = R[kid[row]]; o = int(orb[v]); d[o] = d.get(o, 0) + int(Fr(pr) * SC)
    first = {}
    for v in range(nv): first.setdefault(int(fx.vw[v]), []).append(v)
    for W in keep:
        inner = (W >> 1) & (fx.nw - 1); ol = W >> (n + 1); orr = W & 1
        i0 = (inner >> (n - 1)) & 1; i1 = (inner >> (n - 2)) & 1; j1 = (inner >> 1) & 1; j0 = inner & 1
        c = int(Fr(LT[ol][i0][i1]) * Fr(RT[j1][j0][orr]) * SC); d = R[kid[W]]
        for v in first[inner]:
            o = int(orb[v]); d[o] = d.get(o, 0) - c
    dr = {}                                                      # split flow = merge flow
    for v in range(nv):
        o = int(orb[v]); dr[o] = dr.get(o, 0) + (1 if fx.vname[v][0] == 'S' else -1 if fx.vname[v] == 'M' else 0)
    R.append(dr); R.append({o: int(size[o]) for o in range(no)})   # last row: total flow = 1
    rhs = [0] * (len(R) - 1) + [1]
    nr = len(R); log(f'orbits {no}, rows {nr}')
    Cf = np.zeros((nr, no))
    for k, d in enumerate(R):
        for o, c in d.items(): Cf[k, o] = c
    Cf /= np.abs(Cf).max(axis=1, keepdims=True).clip(min=1)
    log(f'float check of the exact rows at the float optimum: {np.abs(Cf @ xstar - np.array(rhs) / np.abs(np.array([max(1, max(abs(c) for c in d.values())) for d in R]))).max():.2e}')
    Rq, piv = sl.qr(Cf, mode='r', pivoting=True)
    dg = np.abs(np.diag(Rq[0] if isinstance(Rq, tuple) else Rq)); rank = int((dg > 1e-9 * dg[0]).sum())
    B = np.sort(piv[:rank]); log(f'rank {rank}, free {no - rank}')
    Rq2, piv2 = sl.qr(Cf[:, B].T, mode='r', pivoting=True)
    rowsel = np.sort(piv2[:rank])
    Af = Cf[np.ix_(rowsel, B)]; log(f'basis chosen; condition number of the scaled basis matrix {np.linalg.cond(Af):.2e}')
    Bset = set(int(b) for b in B); pos = {int(b): k for k, b in enumerate(B)}
    xN = {o: int(round(xstar[o] * 2 ** FB)) for o in range(no) if o not in Bset}
    assert all(v > 0 for v in xN.values())
    # A y = b with y = 2^FB x_B (integers everywhere)
    Arows = []; b = []
    for k in rowsel:
        d = R[k]; Arows.append([(pos[o], c) for o, c in d.items() if o in Bset])
        b.append(rhs[k] * 2 ** FB - sum(c * xN[o] for o, c in d.items() if o not in Bset))
    A = np.zeros((rank, rank))
    for k, row in enumerate(Arows):
        for j, c in row: A[k, j] = c
    lu = sl.lu_factor(A); log('LU done')
    # iterative refinement with exact residuals: y ~ X / 2^K
    X = [0] * rank; K = 0; res = list(b); kb = 30 if np.linalg.cond(Af) < 1e5 else 20
    target = None; steps = 0
    def solve_float(vec):
        m = max(abs(v) for v in vec)
        if m == 0: return np.zeros(rank), 0
        e = max(0, m.bit_length() - 900)
        return sl.lu_solve(lu, np.array([float(v >> e) if v >= 0 else -float((-v) >> e) for v in vec])), e
    Dcand = None
    while True:
        y, e = solve_float(res)
        # y * 2^e ~ A^-1 res; new digits: round(y 2^e 2^kb)
        inc = [int(round(float(v) * 2 ** 40)) << (e + kb) >> 40 if e + kb >= 40 else int(round(float(v) * 2 ** (e + kb))) for v in y]
        X = [(x << kb) + i for x, i in zip(X, inc)]; K += kb
        res = [(rv << kb) - sum(c * inc[j] for j, c in row) for rv, row in zip(res, Arows)]
        steps += 1
        assert K < 80 * rank + 4000, 'no rational reconstruction: stopping'
        if steps % 250 == 0:
            log(f'  refinement step {steps}: K = {K} bits, residual bits {max(abs(v) for v in res).bit_length()}')
            # try to reconstruct the common denominator from a few components
            D = reconstruct_den(X, K, rank)
            if D is not None:
                N = [round_div(x * D, K) for x in X]
                if all(sum(c * N[j] for j, c in row) == bv * D for row, bv in zip(Arows, b)):
                    log(f'  exact solution found: denominator has {D.bit_length()} bits'); break
    # assemble: flow of orbit o = num[o] / (D 2^FB)
    num = [0] * no
    for o in range(no): num[o] = N[pos[o]] if o in Bset else xN[o] * D
    den = D * 2 ** FB
    assert all(v > 0 for v in num), 'a flow is not positive'
    for d, rv in zip(R, rhs): assert sum(c * num[o] for o, c in d.items()) == rv * den
    log('all orbit rows hold exactly; all flows positive')
    cert = {'r': r, 's': (s.numerator, s.denominator), 't': (t.numerator, t.denominator), 'den': den,
            'flows': {(int(w), nm): num[orb[v]] for v, (w, nm) in enumerate(zip(fx.vw, fx.vname))}}
    with open(out, 'wb') as f: pickle.dump(cert, f, protocol=4)
    hx = -sum(float(Fr(num[orb[v]], den)) * math.log(float(Fr(num[orb[v]], den)) / float(sum(Fr(num[orb[u]], den) for u in first[int(fx.vw[v])]))) for v in range(nv))
    log(f'written {out}: float entropy of the exact flows: exp(2h) = {math.exp(2*hx):.9f}')

def round_div(a, K):
    return (a + (1 << (K - 1))) >> K

def cf_den(num, K, bound):
    """denominator q <= bound of the continued-fraction convergent of num / 2^K that is closest (last one below bound)."""
    a, bq = num, 1 << K; q0, q1 = 1, 0
    while bq:
        t = a // bq; a, bq = bq, a - t * bq
        q0, q1 = q1, q0 + t * q1
        if q1 > bound: return q0
    return q1

def reconstruct_den(X, K, rank):
    """common denominator D of the rationals approximated by X[i] / 2^K (None if K is still too small)."""
    bound = 1 << (K // 2 - 8); D = 1
    for i in range(rank):
        z = X[i] * D; rem = z - (round_div(z, K) << K)
        if abs(rem).bit_length() <= K // 2: continue                 # already an integer up to the refinement error
        q = cf_den(z, K, bound)
        if q <= 1: return None
        D *= q
        if D > bound: return None
        z = X[i] * D; rem = z - (round_div(z, K) << K)
        if abs(rem).bit_length() > K // 2: return None
    return D

if __name__ == '__main__': main()
