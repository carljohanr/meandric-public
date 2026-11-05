"""flow_opt.py (1367): float search for reversible stationary flows on the two-stack machine.

Not part of any proof: it only searches for good parameters.  The certificate is checked exactly by
check_flow_a.py and check_flow_b.py.

Window w = (a, b, c): a, b = the top r symbols of the two stacks (bit k-1 of the integer = k-th symbol
from the top), c = parity.  Reference law beta: a reversible Markov chain of memory m on {0,1}
(m = 0: iid with P(0) = theta; m = 1: P(next = x | prev = y) = T[y][x]).  Free variables
X_A(w) > 0 on windows with a_1 = 1 - c (pop A allowed), X_B likewise.  Detailed balance:
  F(w, popA) = X_A(w),   F(w, pushA) = X_A(pushA(w)) * q(c | a_1..a_m),
  J(w, move) = beta(a) beta(b) / 2 * F(w, move),   h = sum J log(nu/J) / sum J.
Usage: python flow_opt.py r m [theta | t01 t10] [maxiter]
"""
import sys, numpy as np
from scipy.optimize import minimize


def setup(r, m, par):
    mask = (1 << r) - 1
    A = np.arange(1 << r)
    bits = ((A[:, None] >> np.arange(r)[None, :]) & 1)          # bits[a, k] = (k+1)-th symbol from top
    if m == 0:
        th = par[0]
        p0 = np.where(bits == 0, th, 1 - th)
        beta = p0.prod(axis=1)
        def q(c, a):                                             # P(symbol above = c | symbols below)
            return np.where(c == 0, th, 1 - th) * np.ones_like(a, dtype=float)
    else:
        t01, t10 = par                                            # P(0->1), P(1->0); reversible (2 states)
        pi0 = t10 / (t01 + t10); pi = np.array([pi0, 1 - pi0])
        T = np.array([[1 - t01, t01], [t10, 1 - t10]])
        beta = pi[bits[:, 0]].astype(float)
        for k in range(1, r):
            beta = beta * T[bits[:, k - 1], bits[:, k]]
        def q(c, a):                                             # reversible: above-given-below = forward
            return T[a & 1, c]
    return mask, beta, q


def build(r, m, par):
    mask, beta, q = setup(r, m, par)
    N = 1 << (2 * r + 1)
    w = np.arange(N)
    a = w & mask; b = (w >> r) & mask; c = w >> (2 * r)
    bw = beta[a] * beta[b] / 2
    def idx(a_, b_, c_):
        return a_ | (b_ << r) | (c_ << (2 * r))
    popA_ok = (a & 1) != c
    popB_ok = (b & 1) != c
    pushA_t = idx(((a << 1) | c) & mask, b, 1 - c)
    pushB_t = idx(a, ((b << 1) | c) & mask, 1 - c)
    qA = q(c, a); qB = q(c, b)
    return dict(N=N, bw=bw, popA_ok=popA_ok, popB_ok=popB_ok, pushA_t=pushA_t, pushB_t=pushB_t, qA=qA, qB=qB)


def flows(S, XA, XB):
    # XA, XB: arrays of length N (only entries with pop allowed are used)
    bw = S['bw']
    JpA = np.where(S['popA_ok'], bw * XA, 0.0)
    JpB = np.where(S['popB_ok'], bw * XB, 0.0)
    JoA = bw * XA[S['pushA_t']] * S['qA']
    JoB = bw * XB[S['pushB_t']] * S['qB']
    return JoA, JoB, JpA, JpB


def objective(S, y, sel):
    XA = np.zeros(S['N']); XB = np.zeros(S['N'])
    nA = sel[0].sum()
    XA[sel[0]] = np.exp(y[:nA]); XB[sel[1]] = np.exp(y[nA:])
    Js = flows(S, XA, XB)
    nu = sum(Js)
    M = nu.sum()
    H = 0.0; G = []
    for J in Js:
        with np.errstate(divide='ignore', invalid='ignore'):
            lg = np.where(J > 0, np.log(nu / np.where(J > 0, J, 1)), 0.0)
        H += (J * lg).sum()
        G.append(lg)
    phi = H / M
    g = [(lg - phi) / M for lg in G]                       # dphi/dJ for each move
    bw = S['bw']
    gXA = np.where(S['popA_ok'], bw * g[2], 0.0)
    np.add.at(gXA, S['pushA_t'], bw * S['qA'] * g[0])
    gXB = np.where(S['popB_ok'], bw * g[3], 0.0)
    np.add.at(gXB, S['pushB_t'], bw * S['qB'] * g[1])
    grad = np.concatenate([gXA[sel[0]] * XA[sel[0]], gXB[sel[1]] * XB[sel[1]]])
    return -phi, -grad


def optimize(r, m, par, maxiter=3000, y0=None):
    S = build(r, m, par)
    sel = (S['popA_ok'], S['popB_ok'])
    n = sel[0].sum() + sel[1].sum()
    if y0 is None:
        y0 = np.zeros(n)
    res = minimize(lambda y: objective(S, y, sel), y0, jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxiter, maxcor=30, gtol=1e-12, ftol=1e-15))
    return -res.fun, res, S, sel


if __name__ == '__main__':
    r = int(sys.argv[1]); m = int(sys.argv[2])
    par = [float(x) for x in sys.argv[3:3 + (1 if m == 0 else 2)]]
    it = int(sys.argv[3 + len(par)]) if len(sys.argv) > 3 + len(par) else 3000
    h, res, S, sel = optimize(r, m, par, it)
    print('r=%d m=%d par=%s h=%.10f e^{2h}=%.8f iters=%d %s' % (r, m, par, h, np.exp(2 * h), res.nit, res.message))
