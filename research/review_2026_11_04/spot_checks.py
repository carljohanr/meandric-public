"""Small independent checks for a reading review; no certificate replay."""
from collections import defaultdict
from decimal import Decimal, getcontext
from functools import lru_cache
from pathlib import Path
import json

def tape_next(state):
    w, h = state
    if h: yield w, h - 1
    if h + 1 < len(w): yield w, h + 1
    triples = [(1, -1, 1)] if w[h] == 1 else [(1, -1, -1), (-1, -1, 1)]
    for triple in triples:
        yield w[:h] + triple + w[h+1:], h+1
    if w[h] == 1 and 0 < h < len(w)-1:
        s = sum(w[h-1:h+2])
        if s in (-1, 1): yield w[:h-1] + (s,) + w[h+2:], h-1

def tape_count(n):
    root = ((1, -1, 1), 1)
    d = {root: 1}
    for step in range(2*n-2):
        e = defaultdict(int)
        remaining = 2*n-3-step
        for state, count in d.items():
            for nxt in tape_next(state):
                if len(nxt[0])-3 <= 2*remaining:
                    e[nxt] += count
        d = e
    return d.get(root, 0)

@lru_cache(None)
def matchings(points):
    if not points: return ((),)
    result = []
    for j in range(1, len(points), 2):
        for left in matchings(points[1:j]):
            for right in matchings(points[j+1:]):
                result.append(((points[0], points[j]),) + left + right)
    return tuple(result)

def geometric_counts(n):
    total = band = narrow = 0
    for top in matchings(tuple(range(2*n))):
        for bottom in matchings(tuple(range(2*n))):
            adj = [[] for _ in range(2*n)]
            for pairs, sign in [(top, -1), (bottom, 1)]:
                for a, b in pairs:
                    adj[a].append((b, sign)); adj[b].append((a, -sign))
            seen = {0}; todo = [0]
            while todo:
                for b, _ in adj[todo.pop()]:
                    if b not in seen: seen.add(b); todo.append(b)
            if len(seen) != 2*n: continue
            total += 1
            labels = {0: 0}; todo = [0]
            while todo:
                a = todo.pop()
                for b, sign in adj[a]:
                    if b == 2*n-1: continue
                    if b not in labels:
                        labels[b] = labels[a] + sign; todo.append(b)
            assert len(labels) == 2*n-1
            assert labels[2*n-2] == 0
            assert all(abs(labels[i+1]-labels[i]) == 1 for i in range(2*n-2))
            band += min(labels.values()) >= -1 and max(labels.values()) <= 2
            narrow += min(labels.values()) >= -1 and max(labels.values()) <= 1
    return total, band, narrow

def bicubic_count(n):
    root = ((), (), 0)
    d = {root: 1}
    for step in range(2*n):
        e = defaultdict(int)
        remaining = 2*n-1-step
        for (u, l, c), count in d.items():
            nxt = [(u+(c,), l, 1-c), (u, l+(c,), 1-c)]
            if u and u[-1] == 1-c: nxt.append((u[:-1], l, 1-c))
            if l and l[-1] == 1-c: nxt.append((u, l[:-1], 1-c))
            for v in nxt:
                if len(v[0])+len(v[1]) <= remaining: e[v] += count
        d = e
    return d.get(root, 0)

getcontext().prec = 65
lo = Decimal('12.219100857094')
hi = Decimal('12.29693385352921289818')
width = hi-lo
rows = []
for n in range(1, 7):
    geo, band, narrow = geometric_counts(n)
    tape = tape_count(n)
    assert geo == tape
    rows.append(dict(n=n, meanders=tape, band_minus1_2=band,
                     band_minus1_1=narrow, bicubic=bicubic_count(n)))
assert [r['meanders'] for r in rows] == [1,2,8,42,262,1828]
assert [r['band_minus1_2'] for r in rows] == [1,2,7,31,161,932]
assert [r['band_minus1_1'] for r in rows] == [1,2,6,22,92,422]
assert [r['bicubic'] for r in rows] == [2,8,40,228,1424,9520]
x = Decimal('0.92')
for _ in range(20): x -= (4*x**4-2*x-1)/(16*x**3-2)
lam = (1+2*x)*(1+x)**2
alpha = (29+Decimal(145).sqrt())/12
out = dict(scope='Orders 1 through 6 only; no large certificate verification', counts=rows,
    width=str(width), rounded_width=str(Decimal('12.2970')-Decimal('12.2191')),
    reduction_vs_2005_percent=str(100*(1-width/Decimal('1.521'))),
    reduction_vs_Fox_2026_percent=str(100*(1-width/Decimal('1.0239'))),
    factor_vs_2005=str(Decimal('1.521')/width), factor_vs_Fox_2026=str(Decimal('1.0239')/width),
    bicubic_upper_square=str((Decimal(3187011595)/10**9)**2),
    band_growth=str(lam), band_quartic_residual=str(16*lam**4-112*lam**3-520*lam**2-536*lam+25),
    predicted_alpha=str(alpha), continuum_upper_exponent=str(alpha-1+Decimal(5)/7))
Path(__file__).with_name('spot_checks.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
