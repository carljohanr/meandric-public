"""1426 (JB): independent rebuild of the memory-three reference law of Paper B, Section 6.1.

Written from Paper B s6_lower.tex and 1421 SPEC.md Section 3 only. Exact rationals throughout.
Primary reading (Paper B prose): the perturbation minimum is taken over all paths whose middle heights lie in the
orbit of the path's middle heights under reversal and reflection, with paths that leave the strip counting as 0.
Second reading (SPEC Section 3): explicit group labels. Both are built and compared, then compared with law.json.

Usage: python3 law.py <law.json>      (prints a report; exits nonzero on any mismatch)
Importable: build_law() returns a dict with tables keyed by height tuples.
"""
import json
import sys
from fractions import Fraction as Fr

K = 12
UP = [Fr(1), Fr(97, 100), Fr(99, 100), Fr(47, 50), Fr(171, 200), Fr(69, 100), Fr(31, 100), Fr(29, 200),
      Fr(3, 50), Fr(1, 100), Fr(3, 100), Fr(0)]
KAPPA = Fr(4, 5)
ETA = Fr(3, 16)
ZETA = Fr(1, 8)
import os as _os
if _os.environ.get("JB_ZETA"):  # damage test only: a wrong law parameter
    ZETA = Fr(_os.environ["JB_ZETA"])


def paths(k):
    """All height sequences of k heights in {0..K-1} with steps +-1."""
    out = [(h,) for h in range(K)]
    for _ in range(k - 1):
        out = [p + (p[-1] + s,) for p in out for s in (1, -1) if 0 <= p[-1] + s < K]
    return out


def instrip(p):
    return all(0 <= h < K for h in p)


def build_law(up=UP, kappa=KAPPA, eta=ETA, zeta=ZETA, reading="paper"):
    # layer 1+2: edge chain P(a,b,c) = prob of b->c having arrived at b from a
    def P(a, b, c):
        u = up[b]
        if u == 0 or u == 1:
            pu = u
        else:
            o = u / (1 - u) * (kappa if b - a == 1 else 1 / kappa)
            pu = o / (1 + o)
        return pu if c == b + 1 else 1 - pu

    E = paths(2)
    # stationary law of directed edges: solve rho(b,c) = sum_a rho(a,b) P(a,b,c), sum rho = 1, by exact Gaussian
    # elimination (independent of any closed formula); then compare with the detailed-balance formula of SPEC.
    idx = {e: i for i, e in enumerate(E)}
    n = len(E)
    A = [[Fr(0)] * (n + 1) for _ in range(n)]
    for (b, c) in E:
        i = idx[(b, c)]
        A[i][i] += 1
        for a in (b - 1, b + 1):
            if 0 <= a < K:
                A[i][idx[(a, b)]] -= P(a, b, c)
    A[n - 1] = [Fr(1)] * n + [Fr(1)]  # replace one equation by normalization
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    rho = {e: A[idx[e]][n] for e in E}
    # check the full stationarity system (including the equation that was replaced)
    for (b, c) in E:
        assert rho[(b, c)] == sum(rho[(a, b)] * P(a, b, c) for a in (b - 1, b + 1) if 0 <= a < K)
    # SPEC closed form, as a cross-check
    C = [Fr(1)]
    for i in range(1, K - 1):
        C.append(C[-1] * P(i - 1, i, i + 1) / P(i + 1, i, i - 1))
    Z = 2 * sum(C)
    rho_formula = {(a, b): C[min(a, b)] / Z for (a, b) in E}
    assert rho_formula == rho, "edge law: Gaussian solve differs from SPEC closed form"

    pi = {(a, b, c): rho[(a, b)] * P(a, b, c) for (a, b, c) in paths(3)}
    q0 = {(a, b, c, d): pi[(a, b, c)] * P(b, c, d) for (a, b, c, d) in paths(4)}

    def sgn(a, b, d, c):  # +1 if first step up and last step down, or first down and last up
        return (1 if a < b else -1) * (1 if d < c else -1)

    # quartet layer
    q = {}
    for path in paths(4):
        a, b, c, d = path
        if reading == "paper":
            orbit = {(b, c), (c, b), (K - 1 - b, K - 1 - c), (K - 1 - c, K - 1 - b)}
            vals = []
            for (bb, cc) in orbit:
                for aa in (bb - 1, bb + 1):
                    for dd in (cc - 1, cc + 1):
                        pp = (aa, bb, cc, dd)
                        vals.append(q0[pp] if instrip(pp) else Fr(0))
            m = min(vals)
        else:  # SPEC group label
            if 1 <= b <= K - 2 and 1 <= c <= K - 2:
                g = min(b, c, K - 1 - b, K - 1 - c)
                m = min(v for (aa, bb, cc, dd), v in q0.items()
                        if 1 <= bb <= K - 2 and 1 <= cc <= K - 2 and min(bb, cc, K - 1 - bb, K - 1 - cc) == g)
            else:
                m = Fr(0)
        q[path] = q0[path] + eta * m * sgn(a, b, d, c)
    F2 = {(b, c, d, e): q[(b, c, d, e)] / pi[(b, c, d)] for (b, c, d, e) in paths(4)}
    p0 = {pp: q[pp[:4]] * F2[pp[1:]] for pp in paths(5)}
    p = {}
    for path in paths(5):
        a, b, c, d, e = path
        if reading == "paper":
            orbit = {(b, c, d), (d, c, b), (K - 1 - b, K - 1 - c, K - 1 - d), (K - 1 - d, K - 1 - c, K - 1 - b)}
            vals = []
            for (bb, cc, dd) in orbit:
                for aa in (bb - 1, bb + 1):
                    for ee in (dd - 1, dd + 1):
                        pp = (aa, bb, cc, dd, ee)
                        vals.append(p0[pp] if instrip(pp) else Fr(0))
            m = min(vals)
        else:
            if 1 <= b <= K - 2 and 1 <= d <= K - 2:
                def grp(bb, cc, dd):
                    return min((bb, cc, dd), (dd, cc, bb), (K - 1 - bb, K - 1 - cc, K - 1 - dd),
                               (K - 1 - dd, K - 1 - cc, K - 1 - bb))
                G = grp(b, c, d)
                m = min(v for pp, v in p0.items() if 1 <= pp[1] <= K - 2 and 1 <= pp[3] <= K - 2
                        and grp(pp[1], pp[2], pp[3]) == G)
            else:
                m = Fr(0)
        p[path] = p0[path] + zeta * m * sgn(a, b, e, d)
    return dict(rho=rho, pi=pi, q0=q0, q=q, p0=p0, p=p, P=P)


def beta_heights(law, hs):
    """Cylinder probability of a height path hs (len >= 2), exact."""
    k = len(hs) - 1
    if k == 1:
        return law["rho"][tuple(hs)]
    if k == 2:
        return law["pi"][tuple(hs)]
    q, p = law["q"], law["p"]
    v = q[tuple(hs[0:4])]
    for i in range(4, k + 1):
        v = v * p[tuple(hs[i - 4:i + 1])] / q[tuple(hs[i - 4:i])]
    return v


def properties(law):
    rho, pi, q, p = law["rho"], law["pi"], law["q"], law["p"]
    rep = {}
    rep["counts"] = [len(rho), len(pi), len(q), len(p)]
    rep["min_quintet"] = min(p.values())
    rep["all_quintets_positive"] = all(v > 0 for v in p.values())
    ok = True
    for t in pi:
        ok &= sum(q[t + (t[-1] + s,)] for s in (1, -1) if (t + (t[-1] + s,)) in q) == pi[t]
        ok &= sum(q[(t[0] + s,) + t] for s in (1, -1) if ((t[0] + s,) + t) in q) == pi[t]
    for t in q:
        ok &= sum(p[t + (t[-1] + s,)] for s in (1, -1) if (t + (t[-1] + s,)) in p) == q[t]
        ok &= sum(p[(t[0] + s,) + t] for s in (1, -1) if ((t[0] + s,) + t) in p) == q[t]
    for t in rho:
        ok &= sum(pi[t + (t[-1] + s,)] for s in (1, -1) if (t + (t[-1] + s,)) in pi) == rho[t]
        ok &= sum(pi[(t[0] + s,) + t] for s in (1, -1) if ((t[0] + s,) + t) in pi) == rho[t]
    rep["marginals_consistent"] = ok
    rep["total_mass_1"] = sum(rho.values()) == 1
    sym = all(p[t] == p[t[::-1]] and p[t] == p[tuple(K - 1 - h for h in t)] for t in p)
    sym &= all(q[t] == q[t[::-1]] and q[t] == q[tuple(K - 1 - h for h in t)] for t in q)
    rep["reversal_reflection_symmetric"] = sym
    # irreducibility of the chain on states = 4 consecutive heights, transitions with positive quintet prob
    states = list(q)
    succ = {s: [t[1:] for t in p if t[:4] == s and p[t] > 0] for s in states}
    def reach(s0):
        seen = {s0}; st = [s0]
        while st:
            s = st.pop()
            for t in succ[s]:
                if t not in seen:
                    seen.add(t); st.append(t)
        return seen
    rep["irreducible"] = all(len(reach(s)) == len(states) for s in states)
    rep["all_quartets_positive"] = all(v > 0 for v in q.values())
    return rep


def compare_json(law, path):
    d = json.load(open(path))
    mism = []
    assert d["K"] == K and d["memory"] == 3
    if [Fr(x) for x in d["up"]] != UP: mism.append("up")
    if Fr(d["kappa"]) != KAPPA: mism.append("kappa")
    if Fr(d["eta"]) != ETA: mism.append("eta")
    if Fr(d["zeta"]) != ZETA: mism.append("zeta")
    n = 0
    for key, tab in (("rho", law["rho"]), ("pi", law["pi"]), ("quartet", law["q"]), ("quintet", law["p"])):
        ent = {tuple(e[:-1]): Fr(e[-1]) for e in d[key]}
        if set(ent) != set(tab):
            mism.append(key + ": key sets differ")
        for k2, v in tab.items():
            n += 1
            if ent.get(k2) != v:
                mism.append(f"{key}{k2}: file {ent.get(k2)} rebuilt {v}")
    return n, mism


if __name__ == "__main__":
    lawP = build_law(reading="paper")
    lawS = build_law(reading="spec")
    same = all(lawP[k] == lawS[k] for k in ("rho", "pi", "q", "p"))
    print("paper reading == SPEC group reading:", same)
    rep = properties(lawP)
    for k, v in rep.items():
        print(k, float(v) if isinstance(v, Fr) else v)
    n, mism = compare_json(lawP, sys.argv[1])
    print("law.json entries compared:", n, "mismatches:", len(mism))
    for m in mism[:20]:
        print("  ", m)
    # the reading without the off-strip-as-zero rule (minimum over in-strip paths of the orbit only)
    bad = not same or mism or not all(v for k, v in rep.items() if k not in ("counts", "min_quintet"))
    sys.exit(1 if bad else 0)
