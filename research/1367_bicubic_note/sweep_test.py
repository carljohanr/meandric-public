"""sweep_test.py (1367): exact test of the sweep lemma on small states.
(1) The closing strategy of PROOF_LOWER.md (Lemma S) reaches the empty state from every state with
    N0 - N1 = c, in at most 5|U| + 3|L| moves, and every move it makes is a legal move.
(2) Breadth-first search from e over states of size <= n+2 reaches exactly the states with N0 - N1 = c.
(3) The transfer is symmetric: y is a successor of x iff x is a successor of y.
Stacks are strings, bottom first."""
import itertools, sys
from collections import deque

def succ(x):
    U, L, c = x; out = []
    out.append((U + str(c), L, 1 - c)); out.append((U, L + str(c), 1 - c))
    if U and U[-1] != str(c): out.append((U[:-1], L, 1 - c))
    if L and L[-1] != str(c): out.append((U, L[:-1], 1 - c))
    return out

def charge(s): return s.count('0') - s.count('1')

def strategy(x):
    U, L, c = x; moves = 0; path = [x]
    while (U, L, c) != ('', '', 0):
        cs = str(c)
        if U and U[-1] != cs: nxt = [(U[:-1], L, 1 - c)]
        elif L and L[-1] != cs: nxt = [(U, L[:-1], 1 - c)]
        elif any(ch != cs for ch in U):          # transfer top of U to L: push c on L, pop U
            nxt = [(U, L + cs, 1 - c), (U[:-1], L + cs, c)]
        else:                                     # transfer top of L to U
            assert L and L[-1] == cs, (U, L, c)
            nxt = [(U + cs, L, 1 - c), (U + cs, L[:-1], c)]
        for y in nxt:
            assert y in succ(path[-1]), (path[-1], y)
            path.append(y); moves += 1
        U, L, c = path[-1]
    return moves

n = int(sys.argv[1])
worst = 0; cnt = 0
for a in range(n + 1):
    for b in range(n + 1 - a):
        for U in map(''.join, itertools.product('01', repeat=a)):
            for L in map(''.join, itertools.product('01', repeat=b)):
                for c in (0, 1):
                    x = (U, L, c)
                    for y in succ(x): assert x in succ(y)
                    if charge(U) + charge(L) != c: continue
                    mv = strategy(x); cnt += 1
                    assert mv <= 5 * a + 3 * b, (x, mv)
                    if a + b: worst = max(worst, mv / (a + b))
print('strategy: %d states with charge = c and size <= %d closed, worst moves/size %.3f (bound 5)' % (cnt, n, worst))
# BFS
seen = {('', '', 0)}; dq = deque(seen)
while dq:
    x = dq.popleft()
    for y in succ(x):
        if len(y[0]) + len(y[1]) <= n + 2 and y not in seen: seen.add(y); dq.append(y)
bad = [x for x in seen if charge(x[0]) + charge(x[1]) != x[2]]
small = [x for x in seen if len(x[0]) + len(x[1]) <= n]
print('BFS: %d states of size <= %d reached, %d violate the charge rule; %d of size <= %d (strategy set has %d)' % (len(seen), n + 2, len(bad), len(small), n, cnt))
