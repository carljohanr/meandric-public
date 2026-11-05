#!/usr/bin/env python3
"""805: verifier for a stationary-flow certificate (PROOF.md, Theorem F).  Shares no code with flow*.py / exact_cert.py.

Input (pickle): r; s, t (far tape: order-2 chain of the letters read left to right,
P(+ | - -) = s, P(+ | - +) = t, P(+ | + -) = 1 - t, P(+ | + +) = 1 - s); a common denominator `den`; and for every
window w of 2r + 1 letters (integer code, leftmost letter = highest bit, bit 1 = letter +) and every move label m
that is applicable at w a non-negative integer flows[(w, m)]:  J(w, m) = flows[(w, m)] / den.

Checks, all in exact integer / rational arithmetic except the last:
 (V0) the far-tape law: b3 = stationary 3-gram law of the chain (verified), all eight 3-grams positive;
 (V1) J >= 0, J(w, 'L') > 0 for every w, the keys are exactly the applicable (window, move) pairs, total flow 1;
 (V2) stationarity: start from the law  nuW(w) x (far tape continued by the chain), nuW(w) = sum_m J(w, m);
      write the two far letters next to the window on each side explicitly (tape of 2r + 5 letters), make the
      move ON THE TAPE with the machine of 637/641, read the window of radius r + 1 around the new head, and add
      up; the result must equal  nuW(inner window) x P(outer left letter) x P(outer right letter)  for all
      2^(2r+3) windows of radius r + 1;
 (V3) zero drift: sum J * (change of the number of letters left of the head) = 0, same on the right;
 (V4) entropy h = sum_w sum_m J log(nuW(w) / J), lower bound by interval arithmetic (mpmath.iv, 200 bits).
Then (Theorem F)  mu >= exp(2h).

usage: verify_cert.py cert.pkl
"""
import sys, pickle, itertools, time
from fractions import Fraction as Fr
from mpmath import iv

RECORD = Fr(12053547659445446401, 10 ** 18)        # 603: (3471821951 / 10^9)^2

# ---- the tape machine (637 Proposition T; 641 machine.py), on a word (tuple of +1 / -1) with head index h -------
def machine_moves(word, h):
    """list of (label, new word, new head).  The head must have both neighbors."""
    assert 0 < h < len(word) - 1
    x, s, y = word[h - 1], word[h], word[h + 1]
    out = [('L', word, h - 1), ('R', word, h + 1)]
    if s == 1:
        out.append(('S1', word[:h] + (1, -1, 1) + word[h + 1:], h + 1))            # [+] -> + [-] +
    else:
        out.append(('S2', word[:h] + (1, -1, -1) + word[h + 1:], h + 1))           # [-] -> + [-] -
        out.append(('S3', word[:h] + (-1, -1, 1) + word[h + 1:], h + 1))           # [-] -> - [-] +
    if s == 1 and not (x == 1 and y == 1):
        out.append(('M', word[:h - 1] + (x + y + 1,) + word[h + 2:], h - 1))       # x [+] y -> [x + y + 1]
    return out

def check_machine(mmax=7):
    """the machine above, used on finite tapes, counts meanders: closed walks of length 2m - 2 from + [-] +."""
    def fin(word, h):
        w2 = (0,) + word + (0,); out = []
        for m, new, nh in machine_moves(w2, h + 1):
            if m in ('L', 'R') and new[nh] == 0: continue
            if m == 'M' and (w2[h] == 0 or w2[h + 2] == 0): continue
            out.append((new[1:-1], nh - 1))
        return out
    root = ((1, -1, 1), 1); cur = {root: 1}; res = []
    for step in range(2 * mmax - 1):
        if step % 2 == 0: res.append(cur.get(root, 0))
        nxt = {}
        for (word, h), c in cur.items():
            if (len(word) - 3) // 2 > 2 * mmax - 2 - step: continue
            for st in fin(word, h): nxt[st] = nxt.get(st, 0) + c
        cur = nxt
    assert res == [1, 2, 8, 42, 262, 1828, 13820][:mmax], res
    print('(M)  machine check: closed walks from the root =', res)

def code(letters):
    c = 0
    for e in letters: c = (c << 1) | (1 if e == 1 else 0)
    return c

def main():
    t0 = time.time()
    check_machine()
    cert = pickle.load(open(sys.argv[1], 'rb'))
    r = cert['r']; s = Fr(*cert['s']); t = Fr(*cert['t']); den = cert['den']; flows = cert['flows']
    n = 2 * r + 1
    assert r >= 3 and 0 < s < 1 and 0 < t < 1 and den > 0
    # (V0) far tape: transition a[(x, y)] = P(next = + | x, y); stationary pair law b2; 3-gram law b3
    a = {(-1, -1): s, (-1, 1): t, (1, -1): 1 - t, (1, 1): 1 - s}
    u = t / (2 * (s + t)); v = s / (2 * (s + t))
    b2 = {(-1, -1): u, (1, 1): u, (-1, 1): v, (1, -1): v}
    b3 = {(x, y, z): b2[x, y] * (a[x, y] if z == 1 else 1 - a[x, y]) for x in (1, -1) for y in (1, -1) for z in (1, -1)}
    assert sum(b3.values()) == 1 and all(p > 0 for p in b3.values())
    for y in (1, -1):
        for z in (1, -1):                       # stationarity: the law of the last two letters is b2 again
            assert sum(b3[x, y, z] for x in (1, -1)) == b2[y, z]
    def right_of(x, y, z): return b3[x, y, z] / b2[x, y]                 # P(letter z right of the pair x y)
    def left_of(f, x, y): return b3[f, x, y] / (b3[1, x, y] + b3[-1, x, y])   # P(letter f left of the pair x y)
    print(f'(V0) far tape ok: s = {s}, t = {t}, pair law = {u}, {v}')
    # (V1)
    windows = list(itertools.product((1, -1), repeat=n))
    nuW = {}; keys = set(); total = 0
    for w in windows:
        cw = code(w); labels = [m for m, _, _ in machine_moves(w, r)]
        tot = 0
        for m in labels:
            J = flows[(cw, m)]; assert isinstance(J, int) and J >= 0
            if m == 'L': assert J > 0
            tot += J; keys.add((cw, m))
        nuW[w] = tot; total += tot
    assert keys == set(flows.keys()) and total == den
    print(f'(V1) {len(flows)} flows, all >= 0, all L-flows > 0, total flow = 1 exactly; denominator has {den.bit_length()} bits')
    # (V2) one step from the law nuW x chain; everything as integers over the denominator den * Q
    Q = 1
    for p in list(b3.values()) + list(b2.values()): Q = Q * p.denominator // __import__('math').gcd(Q, p.denominator)
    probs = set()
    for x in (1, -1):
        for y in (1, -1):
            for z in (1, -1): probs.add(right_of(x, y, z)); probs.add(left_of(z, x, y))
    Q = 1
    for p in probs: Q = Q * p.denominator // __import__('math').gcd(Q, p.denominator)
    Q4 = Q ** 4
    inflow = {}
    da_sum = 0; db_sum = 0
    for w in windows:
        cw = code(w)
        for f1 in (1, -1):
            p1 = left_of(f1, w[0], w[1])
            for f2 in (1, -1):
                p2 = left_of(f2, f1, w[0])
                for g1 in (1, -1):
                    q1 = right_of(w[-2], w[-1], g1)
                    for g2 in (1, -1):
                        q2 = right_of(w[-1], g1, g2)
                        wt = p1 * p2 * q1 * q2 * Q4; assert wt.denominator == 1; wt = wt.numerator
                        tape = (f2, f1) + w + (g1, g2); h = r + 2
                        for m, new, nh in machine_moves(tape, h):
                            J = flows[(cw, m)]
                            if J == 0: continue
                            assert nh - (r + 1) >= 0 and nh + r + 1 < len(new)
                            Wp = new[nh - (r + 1): nh + r + 2]
                            inflow[Wp] = inflow.get(Wp, 0) + J * wt
        for m, new, nh in machine_moves(w, r):
            J = flows[(cw, m)]
            da_sum += J * (nh - r); db_sum += J * ((len(new) - 1 - nh) - (n - 1 - r))
    bad = 0
    for Wp in itertools.product((1, -1), repeat=n + 2):
        inner = Wp[1:-1]
        want = nuW[inner] * left_of(Wp[0], inner[0], inner[1]) * right_of(inner[-2], inner[-1], Wp[-1]) * Q4
        assert want.denominator == 1
        if inflow.get(Wp, 0) != want.numerator: bad += 1
    assert bad == 0, f'{bad} stationarity equations fail'
    print(f'(V2) stationarity: all {2 ** (n + 2)} equations hold exactly')
    assert da_sum == 0 and db_sum == 0
    print('(V3) zero drift of both side lengths: exact')
    # (V4) entropy
    iv.prec = 200
    D = iv.mpf(den); h = iv.mpf(0)
    for w in windows:
        cw = code(w); nu = iv.mpf(nuW[w])
        for m, _, _ in machine_moves(w, r):
            J = flows[(cw, m)]
            if J > 0: h += (iv.mpf(J) / D) * iv.log(nu / iv.mpf(J))
    from mpmath import mp, nstr
    mp.prec = 200
    low = lambda x: mp.make_mpf(x._mpi_[0])     # lower end point of an interval, as a number
    e2 = iv.exp(2 * h); rec = iv.mpf(RECORD.numerator) / iv.mpf(RECORD.denominator)
    print(f'(V4) entropy: h >= {nstr(low(h), 30)}   (width of the interval {nstr(mp.make_mpf(h._mpi_[1]) - low(h), 3)})')
    print(f'     exp(h)  >= {nstr(low(iv.exp(h)), 30)}')
    print(f'     exp(2h) >= {nstr(low(e2), 30)}')
    print(f'     record of 603: {RECORD.numerator}/10^18;  exp(2h) - record >= {nstr(low(e2 - rec), 12)}')
    # a clean rational statement: the largest k with exp(2h) >= k / 10^9
    k = int(low(e2 * 10 ** 9)); assert low(e2 - iv.mpf(k) / 10 ** 9) >= 0
    print(f'ACCEPTED: mu >= {k}/10^9 = {k / 1e9:.9f}   [{time.time() - t0:.0f}s]')

if __name__ == '__main__': main()
