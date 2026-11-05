"""Referee 1251: own exact checker for the paired reciprocal SOS certificate (written from my own statement,
see STATEMENT.md).  Reads witness JSON as DATA only; imports nothing from 1007/1201/1005.

Words are strings over '0','1' ('1' = step +, '0' = step -).  A state is (s,h).
Weights are kept as integers: every table entry F,E,G is a dyadic rational with denominator dividing 2^36
(asserted on load), stored as numerator over 2^T, T=36.
W(s,h)  = U[window] * prod F[gram]   with U integer, F=Fint/2^T  -> (Wint, nF)
Z(s,h)  = W(complement(s),h)         (complement applied to the padded word when padding is used)
"""
import json
from fractions import Fraction as Q

T = 36
TWO = 1 << T


def load(path):
    return loadd(json.load(open(path)))


def loadd(d):
    c = Cert()
    c.r, c.k, c.q, c.d = d['r'], d['k'], d['q'], d['d']
    c.U = list(d['U'])
    assert all(type(x) is int and x > 0 for x in c.U)

    def toint(x):
        f = Q(x)
        assert f > 0
        n = f * TWO
        assert n.denominator == 1
        return int(n)

    def topair(x):
        f = Q(x)
        assert f > 0
        return (f.numerator, f.denominator)
    c.F = [toint(x) for x in d['F']]      # F must be dyadic with denominator | 2^36 (asserted)
    c.E = [topair(x) for x in d['E']]     # E, G arbitrary positive rationals
    c.G = [topair(x) for x in d['G']]
    c.fl = {int(a): tuple(b) for a, b in d['field_lookup'].items()}
    c.el = {int(a): b for a, b in d['edit_lookup'].items()}
    assert len(c.U) == 1 << (2 * c.r + 1) and len(c.F) == 1 << c.k
    c.raw = d
    return c


class Cert:
    # ---- weights -------------------------------------------------------------------------------------------
    def W(self, s, h, adj=False):
        if adj:
            s = s.translate(_CMP)
        v = self.U[int(s[h - self.r:h + self.r + 1], 2)]
        k = self.k
        n = 0
        F = self.F
        for side in (s[:h], s[h + 1:][::-1]):
            for j in range(len(side) - k + 1):
                v *= F[int(side[j:j + k], 2)]
                n += 1
        return v, n

    # ---- fields: returns (a,b) meaning factor a/b, for the step i -> j ------------------------------------
    def field(self, s, h, t, hj):
        if s == t:  # head shift
            if self.q == 0:
                return 1, 1
            q = self.q
            begin = h - q + (1 if hj > h else 0)
            cut = s[begin:begin + 2 * q]
            assert len(cut) == 2 * q
            idx, sg = self.fl[int(cut, 2)]
            if idx < 0:
                a, b = 1, 1
            elif sg == 1:
                a, b = self.E[idx]
            else:
                b, a = self.E[idx]
            return (a, b) if hj > h else (b, a)
        # edit edge
        if len(s) < len(t):
            short, sh, long_, lh, up = s, h, t, hj, True
        else:
            short, sh, long_, lh, up = t, hj, s, h, False
        d = self.d
        window = short[sh - d:sh + d + 1]
        assert len(window) == 2 * d + 1 and lh == sh + 1
        tri = long_[lh - 1:lh + 2]
        a, b = self.G[self.el[8 * int(window, 2) + int(tri, 2)]]
        return (a, b) if up else (b, a)


_CMP = str.maketrans('01', '10')


# ---- the machine, from my own statement ------------------------------------------------------------------
def out(t):
    s, h = t
    n = len(s)
    if h > 0:
        yield (s, h - 1)
    if h < n - 1:
        yield (s, h + 1)
    if s[h] == '1':
        yield (s[:h] + '101' + s[h + 1:], h + 1)
    else:
        yield (s[:h] + '100' + s[h + 1:], h + 1)
        yield (s[:h] + '001' + s[h + 1:], h + 1)
    if 0 < h < n - 1:
        m = {'010': '0', '011': '1', '110': '1'}.get(s[h - 1:h + 2])
        if m is not None:
            yield (s[:h - 1] + m + s[h + 2:], h - 1)


def inc(t):
    """predecessors j with t in out(j)"""
    s, h = t
    n = len(s)
    if h > 0:
        yield (s, h - 1)
    if h < n - 1:
        yield (s, h + 1)
    # j shorter, t obtained from j by a split: t has head in the middle of 101 / 100 / 001
    if 0 < h < n - 1:
        m = {'101': '1', '100': '0', '001': '0'}.get(s[h - 1:h + 2])
        if m is not None:
            yield (s[:h - 1] + m + s[h + 2:], h - 1)
    # j longer, t obtained from j by a merge: j has middle '+' with triple 010 -> 0 ; 011,110 -> 1
    for tri in (('010',) if s[h] == '0' else ('011', '110')):
        yield (s[:h] + tri + s[h + 1:], h + 1)


def physical(t):
    s, h = t
    return 0 <= h < len(s) and s[0] == s[-1] == '1' and s.count('1') == s.count('0') + 1


def pad(t, p):
    s, h = t
    return '0' * p + s + '0' * p, h + p


def row(c, t, padded=False, keep=None):
    """returns (N, D) with the row value c_t = N/(2D)  hmm: value = (1/2)*(N/D).  Exact integers.
    t is the state in its own coordinates; if padded the weights/fields are evaluated on the r+1-padded words.
    keep: predicate on unpadded neighbour states (restriction); None = all (formal row)."""
    p = c.r + 1 if padded else 0
    nbo = [u for u in out(t) if keep is None or keep(u)]
    nbi = [u for u in inc(t) if keep is None or keep(u)]
    tp = pad(t, p) if p else t
    N, D = 0, 1
    for adj, nbs in ((False, nbo), (True, nbi)):
        wi, ni = c.W(tp[0], tp[1], adj)
        for u in nbs:
            up = pad(u, p) if p else u
            wj, nj = c.W(up[0], up[1], adj)
            a, b = c.field(tp[0], tp[1], up[0], up[1])
            # term = wj * 2^(-T nj) / (wi * 2^(-T ni)) * a/b
            num = wj * a
            den = wi * b
            e = ni - nj
            if e >= 0:
                num <<= T * e
            else:
                den <<= T * (-e)
            N = N * den + num * D
            D = D * den
    return N, D


def leq(N, D, Rn, Rd):
    """(1/2) N/D <= Rn/Rd"""
    return N * Rd <= 2 * Rn * D
