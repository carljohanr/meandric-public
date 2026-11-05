#!/usr/bin/env python3
"""805: stationary-flow program (float; data only, NOT a bound).

Process: the tape seen from the head, bulk machine (804 Definition 1), optionally confined to a strip of K levels.
Law nu = (law nuW of the window of radius r around the head) x (far tape = stationary nearest-neighbor Markov path
continued from the two ends of the window; without a strip: independent letters, + with probability alpha).
Kernel: the move depends on the window only.  Unknowns: the flows J(w, m) = nuW(w) p(m | w) >= 0.
For a FIXED far-tape law the conditions are LINEAR in J:
  (S) stationarity: for every window W' of radius r + 1,
        sum over (w, m, fresh letters) that produce W' of J(w, m) P(fresh) = nuW(inner W') P(outer letters of W'),
  (D) zero drift: flow(L) = flow(R), flow(splits) = flow(merges),      (N) total flow 1.
Entropy of the moves given the past (804 Lemma 3 (e)): h = sum_w [ -sum_m J log J + nuW log nuW ], concave in J.
So for a fixed far-tape law the best process of this shape is a concave maximization over a polytope.
Theorem L (804) then gives  growth >= exp(2h)  when the process is ergodic (PROOF.md, Lemma G).

usage: flow.py K r [up-probabilities at levels 1..K-2 | alpha]      (K = 0: no strip)
"""
import sys, math, itertools
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import qr

def strip_ok(lev0, letters, K):
    if not K: return True
    v = lev0
    if not (0 <= v < K): return False
    for c in letters:
        v += c
        if not (0 <= v < K): return False
    return True

def all_windows(n, K):
    res = []
    for lev0 in (range(K) if K else [0]):
        for letters in itertools.product((1, -1), repeat=n):
            if strip_ok(lev0, letters, K): res.append((lev0, letters))
    return res

class Model:
    def __init__(self, K, r, up):
        """up: for a strip, dict level -> probability that the path steps up from that level (levels 0..K-1,
        up[0] = 1, up[K-1] = 0); without a strip a single number alpha = P(letter +)."""
        self.K, self.r, self.up = K, r, up
        self.wins = all_windows(2 * r + 1, K)
        self.widx = {w: i for i, w in enumerate(self.wins)}
        self.rows = all_windows(2 * r + 3, K)
        self.ridx = {w: i for i, w in enumerate(self.rows)}
        self.vars = []          # (window index, move name)
        self.trans = []         # for each var: list of (row index, probability)
        self.dab = []           # (delta a, delta b)
        for wi, (lev0, e) in enumerate(self.wins):
            for name, d, outs in self.moves(lev0, e):
                self.vars.append((wi, name)); self.dab.append(d)
                if not K: outs = [((0, W[1]), pr) for W, pr in outs]
                self.trans.append([(self.ridx[W], pr) for W, pr in outs if pr > 0])
        self.nv = len(self.vars)

    # fresh letters.  Left: the letter that ENDS at the vertex of level lev; right: the letter that STARTS there.
    def left_fresh(self, lev):
        if not self.K: return [(1, self.up), (-1, 1 - self.up)]
        # reversible birth-death chain: P(previous vertex = lev - 1) = P(step down from lev) = 1 - up[lev]
        return [(1, 1 - self.up[lev]), (-1, self.up[lev])]
    def right_fresh(self, lev):
        if not self.K: return [(1, self.up), (-1, 1 - self.up)]
        return [(1, self.up[lev]), (-1, 1 - self.up[lev])]

    def moves(self, lev0, e):
        r, K = self.r, self.K
        A, s, B = e[:r], e[r], e[r + 1:]
        endlev = lev0 + sum(e)
        res = []
        # L: two fresh letters on the left
        outs = []
        for f1, p1 in self.left_fresh(lev0):
            if p1 == 0: continue
            for f2, p2 in self.left_fresh(lev0 - f1):
                outs.append(((lev0 - f1 - f2, (f2, f1) + e), p1 * p2))
        res.append(('L', (-1, 1), outs))
        outs = []
        for g1, p1 in self.right_fresh(endlev):
            if p1 == 0: continue
            for g2, p2 in self.right_fresh(endlev + g1):
                outs.append(((lev0, e + (g1, g2)), p1 * p2))
        res.append(('R', (1, -1), outs))
        splits = [('S1', (1, -1, 1))] if s == 1 else [('S2', (1, -1, -1)), ('S3', (-1, -1, 1))]
        for name, trip in splits:
            new = A + trip + B
            if strip_ok(lev0, new, K): res.append((name, (1, 1), [((lev0, new), 1.0)]))
        if s == 1 and not (A[-1] == 1 and B[0] == 1):
            z = A[-1] + B[0] + 1
            core = A[:-1] + (z,) + B[1:]
            outs = []
            for f1, p1 in self.left_fresh(lev0):
                if p1 == 0: continue
                for f2, p2 in self.left_fresh(lev0 - f1):
                    for g1, q1 in self.right_fresh(endlev):
                        if q1 == 0: continue
                        for g2, q2 in self.right_fresh(endlev + g1):
                            outs.append(((lev0 - f1 - f2, (f2, f1) + core + (g1, g2)), p1 * p2 * q1 * q2))
            res.append(('M', (-1, -1), outs))
        return res

    def constraints(self):
        """matrix C and vector d with C J = d: (S) for every row, (D), (N)."""
        nr = len(self.rows); C = np.zeros((nr + 3, self.nv)); d = np.zeros(nr + 3)
        first = {}
        for v, (wi, name) in enumerate(self.vars): first.setdefault(wi, []).append(v)
        for v, tr in enumerate(self.trans):
            for ri, pr in tr: C[ri, v] += pr
        for ri, (lev0, W) in enumerate(self.rows):
            inner = (lev0 + W[0] if self.K else 0, W[1:-1])
            pl = dict(self.left_fresh(inner[0]))[W[0]]
            pr = dict(self.right_fresh(inner[0] + sum(inner[1])))[W[-1]]
            for v in first[self.widx[inner]]: C[ri, v] -= pl * pr
        for v, (da, db) in enumerate(self.dab):
            C[nr, v] = da; C[nr + 1, v] = db
        C[nr + 2, :] = 1; d[nr + 2] = 1
        return C, d

    def entropy(self, J):
        h = 0.0; nuW = np.zeros(len(self.wins))
        for v, (wi, _) in enumerate(self.vars): nuW[wi] += J[v]
        for v, (wi, _) in enumerate(self.vars):
            if J[v] > 0: h -= J[v] * math.log(J[v] / nuW[wi])
        return h, nuW

def solve(model, verbose=False, x0=None):
    C, d = model.constraints()
    # independent rows
    Q, Rm, piv = qr(C.T, mode='economic', pivoting=True)
    diag = np.abs(np.diag(Rm)); rank = int((diag > 1e-10 * diag[0]).sum())
    rows = np.sort(piv[:rank]); Cr, dr = C[rows], d[rows]
    nv = model.nv
    win = np.array([wi for wi, _ in model.vars]); nw = len(model.wins)
    def negH(J):
        J = np.maximum(J, 1e-300); nu = np.bincount(win, weights=J, minlength=nw)
        return float(np.sum(J * (np.log(J) - np.log(nu[win]))))
    def grad(J):
        J = np.maximum(J, 1e-300); nu = np.bincount(win, weights=J, minlength=nw)
        return np.log(J) - np.log(nu[win])
    if x0 is None:
        # least-norm solution pushed to the interior
        x0 = np.linalg.lstsq(Cr, dr, rcond=None)[0]; x0 = np.maximum(x0, 1e-4); x0 /= x0.sum()
    res = minimize(negH, x0, jac=grad, method='SLSQP', bounds=[(1e-13, 1)] * nv,
                   constraints=[{'type': 'eq', 'fun': lambda J: Cr @ J - dr, 'jac': lambda J: Cr}],
                   options={'maxiter': 3000, 'ftol': 1e-15})
    J = res.x; h, nuW = model.entropy(J)
    viol = float(np.abs(C @ J - d).max())
    if verbose: print(f'  vars {nv}, rows {len(model.rows)}, rank {rank}, free dim {nv - rank}, iterations {res.nit}, '
                      f'status {res.status}, max violation {viol:.2e}')
    return h, J, nuW, viol, nv - rank

def sym_up(K, params):
    """mirror-symmetric up-probabilities for a strip of K levels from the free parameters (levels 1..floor((K-1)/2))."""
    up = {0: 1.0, K - 1: 0.0}
    inner = list(range(1, K - 1)); half = (len(inner) + 1) // 2
    for j, lev in enumerate(inner[:half]):
        if len(inner) % 2 == 1 and j == half - 1: up[lev] = 0.5
        else: up[lev] = params[j]; up[K - 1 - lev] = 1 - params[j]
    return up

if __name__ == '__main__':
    K = int(sys.argv[1]); r = int(sys.argv[2]); par = [float(x) for x in sys.argv[3:]]
    if K: up = sym_up(K, par); print('up-probabilities:', up)
    else: up = par[0] if par else 0.5; print('alpha =', up)
    m = Model(K, r, up)
    h, J, nuW, viol, free = solve(m, verbose=True)
    print(f'K={K} r={r}: h = {h:.10f}, exp(h) = {math.exp(h):.10f}, exp(2h) = {math.exp(2*h):.10f}')
    flows = {}
    for v, (wi, name) in enumerate(m.vars): flows[name[0]] = flows.get(name[0], 0) + J[v]
    print('  flows by move kind:', {k: round(v, 6) for k, v in flows.items()}, ' min window prob', nuW.min())
