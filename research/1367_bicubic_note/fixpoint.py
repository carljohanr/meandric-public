"""fixpoint.py (1367): fixed-point iteration for the optimal reversible flow (float search only).
At an optimum, log X(v) = weighted mean over the edges e that use X(v) of log(nu(w_e)/coef_e) + const.
Usage: python fixpoint.py r m s iters [save.npz]"""
import sys, os, numpy as np, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flow_opt import build, flows

def entropy(S, XA, XB):
    Js = flows(S, XA, XB); nu = sum(Js); M = nu.sum(); H = 0.0
    for J in Js:
        pos = J > 0
        H += (J[pos] * np.log(nu[pos] / J[pos])).sum()
    return H / M, Js, nu

def run(r, m, par, iters, omega=1.0, XA=None, XB=None, report=200):
    S = build(r, m, par); N = S['N']; bw = S['bw']
    selA, selB = S['popA_ok'], S['popB_ok']
    if XA is None:
        XA = np.where(selA, 1.0, 0.0); XB = np.where(selB, 1.0, 0.0)
    # coefficients of the edges using X_A(v): (v, popA) coef bw[v]; (w, pushA) with pushA_t[w] = v coef bw[w] qA[w]
    cA_pop = np.where(selA, bw, 0.0); cA_push = bw * S['qA']
    cB_pop = np.where(selB, bw, 0.0); cB_push = bw * S['qB']
    totA = cA_pop.copy(); np.add.at(totA, S['pushA_t'], cA_push)
    totB = cB_pop.copy(); np.add.at(totB, S['pushB_t'], cB_push)
    hist = []
    for it in range(iters + 1):
        if it % report == 0 or it == iters:
            h, Js, nu = entropy(S, XA, XB)
            hist.append((it, h)); print('it %6d  h %.12f  e^2h %.10f' % (it, h, np.exp(2 * h)), flush=True)
            if it == iters: break
        Js = flows(S, XA, XB); nu = sum(Js)
        lnu = np.log(nu)
        # weighted mean of log(nu/coef) = log nu - log coef (coef = bw or bw q)
        sA = np.where(selA, cA_pop * (lnu - np.log(np.where(selA, cA_pop, 1))), 0.0)
        np.add.at(sA, S['pushA_t'], cA_push * (lnu - np.log(cA_push)))
        sB = np.where(selB, cB_pop * (lnu - np.log(np.where(selB, cB_pop, 1))), 0.0)
        np.add.at(sB, S['pushB_t'], cB_push * (lnu - np.log(cB_push)))
        nA = np.where(selA, sA / np.where(selA, totA, 1), 0.0)
        nB = np.where(selB, sB / np.where(selB, totB, 1), 0.0)
        lA = np.where(selA, np.log(np.where(selA, XA, 1)), 0.0); lB = np.where(selB, np.log(np.where(selB, XB, 1)), 0.0)
        lA = (1 - omega) * lA + omega * nA; lB = (1 - omega) * lB + omega * nB
        sh = max(lA[selA].max(), lB[selB].max())
        XA = np.where(selA, np.exp(lA - sh), 0.0); XB = np.where(selB, np.exp(lB - sh), 0.0)
    return XA, XB, S, hist

if __name__ == '__main__':
    r, m = int(sys.argv[1]), int(sys.argv[2]); s = float(sys.argv[3]); iters = int(sys.argv[4])
    par = [s] if m == 0 else [s, s]
    om = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
    t0 = time.time()
    XA, XB, S, hist = run(r, m, par, iters, omega=om, report=max(1, iters // 10))
    print('time %.1f s' % (time.time() - t0))
    if len(sys.argv) > 5 and sys.argv[5] != '-':
        np.savez(sys.argv[5], XA=XA, XB=XB, r=r, m=m, s=s)
