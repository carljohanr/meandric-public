#!/usr/bin/env python3
"""805: stationary-flow program for ALL meanders, scalable version (float; data only, NOT a bound).

Same program as flow.py with K = 0 (no strip), far tape = stationary two-state Markov chain of the letters read left
to right: P(next = + | current = +) = pp, P(next = + | current = -) = pm (independent letters: pp = pm = alpha).
The two-state chain is reversible, so a fresh letter on the left has the same conditional law given its right
neighbor.  Windows are bit codes (bit 1 = letter +, leftmost letter = highest bit).
Solver: infeasible-start damped Newton on  min sum J log(J / nuW)  subject to  C J = d  (sparse KKT systems).

usage: flow2.py r [pp pm]        (default pp = pm = 0.5)
"""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'): os.environ[v] = '1'
import sys, math, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

class Flow:
    def __init__(self, r, pp=0.5, pm=0.5, b3=None):
        """b3: optional stationary 3-gram law of the far tape, b3[x][y][z] (bits, left to right); then the far tape is
        the order-2 chain with these 3-grams (needs r >= 2).  Otherwise the order-1 chain (pp, pm)."""
        self.r = r; n = self.n = 2 * r + 1; self.pp, self.pm = pp, pm
        self.nw = 1 << n; self.nrow = 1 << (n + 2)
        if b3 is None:
            T1 = [[1 - pm, pm], [1 - pp, pp]]; st1 = pm / (1 - pp + pm); st = [1 - st1, st1]
            b3 = [[[st[x] * T1[x][y] * T1[y][z] for z in (0, 1)] for y in (0, 1)] for x in (0, 1)]
        b3 = np.array(b3, dtype=object if not isinstance(b3[0][0][0], float) else float)
        self.b3 = b3
        b2 = [[b3[x][y][0] + b3[x][y][1] for y in (0, 1)] for x in (0, 1)]
        # RT[x][y][z] = P(letter z right of the pair (x, y)),  LT[f][x][y] = P(letter f left of the pair (x, y))
        RT = [[[b3[x][y][z] / b2[x][y] for z in (0, 1)] for y in (0, 1)] for x in (0, 1)]
        LT = [[[b3[f][x][y] / (b3[0][x][y] + b3[1][x][y]) for y in (0, 1)] for x in (0, 1)] for f in (0, 1)]
        self.RT, self.LT = RT, LT
        vw, vname, vda, vdb = [], [], [], []
        ri, ci, val = [], [], []
        def bit(code, i, length): return (code >> (length - 1 - i)) & 1          # letter i from the left
        for w in range(self.nw):
            e = [bit(w, i, n) for i in range(n)]
            s = e[r]; x = e[r - 1]; y = e[r + 1]
            def add(name, da, db, outs):
                v = len(vw); vw.append(w); vname.append(name); vda.append(da); vdb.append(db)
                for row, pr in outs:
                    if pr > 0: ri.append(row); ci.append(v); val.append(pr)
            # L: fresh f1 left of e[0], f2 left of f1 (reversible chain: P(f | right neighbor g) = T[g][f])
            outs = []
            for f1 in (0, 1):
                for f2 in (0, 1):
                    outs.append(((f2 << (n + 1)) | (f1 << n) | w, LT[f1][e[0]][e[1]] * LT[f2][f1][e[0]]))
            add('L', -1, 1, outs)
            outs = []
            for g1 in (0, 1):
                for g2 in (0, 1):
                    outs.append(((w << 2) | (g1 << 1) | g2, RT[e[-2]][e[-1]][g1] * RT[e[-1]][g1][g2]))
            add('R', 1, -1, outs)
            left = w >> (r + 1); right = w & ((1 << r) - 1)
            for name, trip in ([('S1', 0b101)] if s == 1 else [('S2', 0b100), ('S3', 0b001)]):
                add(name, 1, 1, [((left << (r + 3)) | (trip << r) | right, 1.0)])
            if s == 1 and not (x == 1 and y == 1):
                z = 1 if (x + y) == 1 else 0                       # (-,-) -> -, mixed -> +
                core = ((left >> 1) << r) | (z << (r - 1)) | (right & ((1 << (r - 1)) - 1))     # 2r - 1 letters
                cw = [bit(core, i, n - 2) for i in range(n - 2)]
                if len(cw) == 1: cw = [cw[0], cw[0]]        # r = 1: only for order-1 chains (second letter unused)
                outs = []
                for f1 in (0, 1):
                    for f2 in (0, 1):
                        for g1 in (0, 1):
                            for g2 in (0, 1):
                                row = (f2 << (n + 1)) | (f1 << n) | (core << 2) | (g1 << 1) | g2
                                outs.append((row, LT[f1][cw[0]][cw[1]] * LT[f2][f1][cw[0]]
                                             * RT[cw[-2]][cw[-1]][g1] * RT[cw[-1]][g1][g2]))
                add('M', -1, -1, outs)
        self.nv = len(vw); self.vw = np.array(vw); self.vname = vname
        self.tri = (ri, ci, val)                      # exact transition entries when b3 is given in fractions
        self.vda = np.array(vda, float); self.vdb = np.array(vdb, float)
        A = sp.csr_matrix((np.array(val, dtype=float), (ri, ci)), shape=(self.nrow, self.nv))
        # right side of (S): nuW(inner W') P(outer left | inner first) P(outer right | inner last)
        rows = np.arange(self.nrow)
        inner = (rows >> 1) & (self.nw - 1)
        ol = rows >> (n + 1); orr = rows & 1
        i0 = (inner >> (n - 1)) & 1; i1 = (inner >> (n - 2)) & 1; j1 = (inner >> 1) & 1; j0 = inner & 1
        LTn = np.array(LT, dtype=float); RTn = np.array(RT, dtype=float)
        coef = LTn[ol, i0, i1] * RTn[j1, j0, orr]
        B = sp.csr_matrix((coef, (rows, inner)), shape=(self.nrow, self.nw))
        W = sp.csr_matrix((np.ones(self.nv), (self.vw, np.arange(self.nv))), shape=(self.nw, self.nv))
        S = A - B @ W
        D = sp.csr_matrix(np.vstack([self.vda, self.vdb, np.ones(self.nv)]))
        self.C = sp.vstack([S, D]).tocsr(); self.d = np.zeros(self.nrow + 3); self.d[-1] = 1.0
        self.Wm = W

    def entropy(self, J):
        nu = self.Wm @ J
        return -float(np.sum(J * (np.log(J) - np.log(nu[self.vw])))), nu

    def solve(self, iters=200, tol=1e-12, verbose=True, J0=None):
        C, d, nv, nc = self.C, self.d, self.nv, self.C.shape[0]
        J = np.full(nv, 1.0 / nv) if J0 is None else J0.copy()
        eps, delta = 1e-10, 1e-9
        t0 = time.time()
        for it in range(iters):
            nu = self.Wm @ J
            g = np.log(J) - np.log(nu[self.vw])
            # Hessian: diag(1/J) - W^T diag(1/nu) W
            H = sp.diags(1.0 / J + eps) - self.Wm.T @ sp.diags(1.0 / nu) @ self.Wm
            res = d - C @ J
            KKT = sp.bmat([[H, C.T], [C, -delta * sp.eye(nc)]], format='csc')
            rhs = np.concatenate([-g, res])
            lu = spl.splu(KKT, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
            sol = lu.solve(rhs)
            sol += lu.solve(rhs - KKT @ sol)
            dJ = sol[:nv]
            neg = dJ < 0
            tmax = np.min(-J[neg] / dJ[neg]) if neg.any() else np.inf
            t = min(1.0, 0.9 * tmax)
            J = J + t * dJ
            h = -float(np.sum(J * g)); rn = float(np.abs(d - C @ J).max()); step = float(np.abs(dJ).max())
            if verbose and (it < 5 or it % 10 == 0):
                print(f'    it {it}: h~{h:.10f} res {rn:.2e} step {step:.2e} t {t:.3f}  [{time.time()-t0:.0f}s]', flush=True)
            if step < 1e-11 and t == 1.0: break
        # remove the residual left by the regularization: least-change projections onto C J = d
        P = spl.splu(sp.bmat([[sp.eye(nv), C.T], [C, -delta * sp.eye(nc)]], format='csc'),
                     permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
        for _ in range(6):
            res = d - C @ J
            J = J + P.solve(np.concatenate([np.zeros(nv), res]))[:nv]
        rn = float(np.abs(d - C @ J).max())
        h, nu = self.entropy(J)
        return h, J, nu, rn, it

if __name__ == '__main__':
    r = int(sys.argv[1]); pp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    pm = float(sys.argv[3]) if len(sys.argv) > 3 else pp
    t0 = time.time(); f = Flow(r, pp, pm)
    print(f'r={r} pp={pp} pm={pm}: windows {f.nw}, flows {f.nv}, rows {f.nrow}  (build {time.time()-t0:.1f}s)', flush=True)
    h, J, nu, rn, it = f.solve()
    kinds = {}
    for v, nm in enumerate(f.vname): kinds[nm[0]] = kinds.get(nm[0], 0.0) + J[v]
    print(f'r={r}: h = {h:.12f}  exp(h) = {math.exp(h):.10f}  exp(2h) = {math.exp(2*h):.10f}   residual {rn:.1e}, '
          f'{it+1} Newton steps, {time.time()-t0:.0f}s')
    print('  flows by kind:', {k: round(float(v), 6) for k, v in kinds.items()}, ' min J', float(J.min()), ' min nuW', float(nu.min()))
    np.save(f'J_r{r}_pp{pp}_pm{pm}.npy', J)
