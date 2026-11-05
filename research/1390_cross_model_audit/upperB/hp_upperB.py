"""1390 (HP): Paper B upper bound, independent check (Lemma 6, Corollary 7, Lemma 8, Proposition 9).

Two implementations written from Paper B Sections 2-4 only (no project checker read):
  * hp_rowsB.c: float filter over all 2^(2r+3) formal rows (moves made literally on padded words);
  * this file: an exact Fraction evaluator, written separately, used (a) on every row the filter flags,
    (b) on random rows to cross-check the two implementations, (c) on actual tapes of the machine, with
    full-tape weights and only the moves the machine really has (a direct test of the reduction).
The JSON encoding is not documented in the project; the reading used here is chosen by option flags and
justified by the paper's symmetries (Section 4.2).  By Proposition 9 ANY positive tables with reciprocal
factors give a valid bound, so the reading can only change the number, not the validity.

Usage:
  python hp_upperB.py prep cand.json tables.bin [fbits=msb|lsb] [ebits=msb|lsb] [gbits=msb|lsb]
  python hp_upperB.py exact cand.json flagged.txt Rnum/Rden [opts]      exact re-evaluation of flagged rows
  python hp_upperB.py sample cand.json N [opts]                          exact vs float on random rows (prints)
  python hp_upperB.py tapes cand.json N maxlen [opts]                    actual tapes vs formal rows
"""
import sys, json, struct, random
from fractions import Fraction as Fr

def load(path, opts):
    d = json.load(open(path))
    r, k, q, dd = d['r'], d['k'], d['q'], d['d']
    fbits = opts.get('fbits', 'msb'); ebits = opts.get('ebits', 'msb'); gbits = opts.get('gbits', 'msb')
    def rev(x, n):
        return int(format(x, '0%db' % n)[::-1], 2)
    nU = 1 << (2 * r + 1)
    Uj = d['U']; assert len(Uj) == nU and all(isinstance(u, int) and u > 0 for u in Uj), 'U not positive integers'
    Fj = [Fr(s) for s in d['F']]; assert len(Fj) == 1 << k and all(x > 0 for x in Fj), 'F not positive'
    Ej = [Fr(s) for s in d['E']]; assert all(x > 0 for x in Ej), 'E not positive'
    Gj = [Fr(s) for s in d['G']]; assert all(x > 0 for x in Gj), 'G not positive'
    fl = d['field_lookup']; el = d['edit_lookup']
    assert len(fl) == 1 << (2 * q), 'field_lookup size'
    U = Uj  # reversal-symmetric (checked below), so the bit order of a window does not matter
    F = Fj if fbits == 'msb' else [Fj[rev(c, k)] for c in range(1 << k)]
    def Eval(cj):
        idx, sg = fl[str(cj)]
        if sg == 0:
            return Fr(1)
        return Ej[idx] if sg == 1 else 1 / Ej[idx]
    E = [Eval(c if ebits == 'msb' else rev(c, 2 * q)) for c in range(1 << (2 * q))]
    G = {}
    nw = 2 * dd + 1
    for key, idx in el.items():
        key = int(key); w, t = key >> 3, key & 7
        if gbits == 'lsb':
            w, t = rev(w, nw), rev(t, 3)
        G[(w << 3) | t] = Gj[idx]
    T = dict(r=r, k=k, q=q, d=dd, U=U, F=F, E=E, G=G, raw=d)
    dmg = opts.get('damage')
    if dmg:
        apply_damage(T, dmg)
    return T

def apply_damage(T, dmg):
    # damage tables at entries used by the stated argmax row 697690 (internal codes)
    r, k, q, dd = T['r'], T['k'], T['q'], T['d']
    w, h = formal_word(T, 697690)
    code = lambda a: int(''.join(map(str, a)), 2)
    if dmg == 'U-':
        T['U'] = list(T['U']); i = code(w[h - r:h + r + 1]); T['U'][i] = T['U'][i] * 999 // 1000
    elif dmg in ('E+', 'E-'):
        E = list(T['E']); i = code(w[h - q + 1:h + q + 1]); M = (1 << (2 * q)) - 1
        f = Fr(1001, 1000) if dmg == 'E+' else Fr(999, 1000)
        rv = lambda x: int(format(x, '0%db' % (2 * q))[::-1], 2)
        for j in {i, i ^ M}:            # keep the stored symmetries: complement-even, reversal-inverse
            E[j] = E[j] * f
        for j in {rv(i), rv(i) ^ M}:
            E[j] = E[j] / f
        T['E'] = E
    elif dmg in ('F+', 'F-'):
        F = list(T['F']); i = code(w[h - k:h])  # gram ending next to the head on the left side
        F[i] = F[i] * (Fr(1001, 1000) if dmg == 'F+' else Fr(999, 1000)); T['F'] = F
    elif dmg in ('G+', 'G-'):
        G = dict(T['G']); wv = code(w[h - dd:h + dd + 1]); tr = (1, 0, 1) if w[h] == 1 else (1, 0, 0)
        key = (wv << 3) | (tr[0] * 4 + tr[1] * 2 + tr[2]); G[key] = G[key] * (Fr(1001, 1000) if dmg == 'G+' else Fr(999, 1000)); T['G'] = G
    elif dmg == 'Fzero':
        F = list(T['F']); F[5] = Fr(0); T['F'] = F
        assert all(x > 0 for x in T['F']), 'F not positive (damaged)'
    else:
        raise SystemExit('unknown damage')
    print('damage applied:', dmg)

def symmetry_report(T):
    r, k, q, dd = T['r'], T['k'], T['q'], T['d']
    def rev(x, n):
        return int(format(x, '0%db' % n)[::-1], 2)
    U = T['U']; n = 2 * r + 1
    badU = sum(1 for i in range(1 << n) if U[i] != U[rev(i, n)])
    E = T['E']; m = 2 * q; M = (1 << m) - 1
    badEc = sum(1 for c in range(1 << m) if E[c] != E[c ^ M])
    badEr = sum(1 for c in range(1 << m) if E[c] * E[rev(c, m)] != 1)
    G = T['G']; nw = 2 * dd + 1; W = (1 << nw) - 1
    badGc = sum(1 for key, g in G.items() if G.get((((key >> 3) ^ W) << 3) | ((key & 7) ^ 7)) != g)
    badGr = sum(1 for key, g in G.items() if G.get((rev(key >> 3, nw) << 3) | rev(key & 7, 3)) != g)
    print(f'symmetry mismatches: U reversal {badU}; E complement {badEc}; E reversal-inverse {badEr}; '
          f'G complement {badGc}; G reversal {badGr}; G entries {len(G)}')

# ---------------- exact, literal evaluation (separately written from the C filter) ----------------
class Machine:
    def __init__(self, T):
        self.__dict__.update(T)
    def W(self, word, h):
        """W of a padded word (list of 0/1, + = 1) with head h, as an exact Fraction."""
        r, k = self.r, self.k
        win = word[h - r:h + r + 1]
        val = Fr(self.U[int(''.join(map(str, win)), 2)])
        left = word[:h]
        right = word[h + 1:][::-1]   # read from the far end toward the head
        for side in (left, right):
            for j in range(len(side) - k + 1):
                val *= self.F[int(''.join(map(str, side[j:j + k])), 2)]
        return val
    def Z(self, word, h):
        return self.W([1 - x for x in word], h)
    def Ecut(self, word, start):
        return self.E[int(''.join(map(str, word[start:start + 2 * self.q])), 2)]
    def Gv(self, short, sh, triple):
        dd = self.d
        w = int(''.join(map(str, short[sh - dd:sh + dd + 1])), 2)
        t = triple[0] * 4 + triple[1] * 2 + triple[2]
        return self.G[(w << 3) | t]
    def moves(self, word, h, real_lo=None, real_hi=None):
        """successors as (theta, word', h'). real_lo/real_hi: indices of the first/last real letter
        (for actual tapes, head moves and merges must stay inside the real tape)."""
        lo = 0 if real_lo is None else real_lo
        hi = len(word) - 1 if real_hi is None else real_hi
        out = []
        if h - 1 >= lo:
            out.append((1 / self.Ecut(word, h - self.q), word, h - 1))
        if h + 1 <= hi:
            out.append((self.Ecut(word, h - self.q + 1), word, h + 1))
        s = word[h]
        trips = [(1, 0, 1)] if s == 1 else [(1, 0, 0), (0, 0, 1)]
        for tr in trips:
            w2 = word[:h] + list(tr) + word[h + 1:]
            out.append((self.Gv(word, h, tr), w2, h + 1))
        if h - 1 >= lo and h + 1 <= hi and s == 1 and not (word[h - 1] == 1 and word[h + 1] == 1):
            x, y = word[h - 1], word[h + 1]
            m = 1 if (2 * x - 1) + 1 + (2 * y - 1) > 0 else 0
            w2 = word[:h - 1] + [m] + word[h + 2:]
            out.append((1 / self.Gv(w2, h - 1, (x, 1, y)), w2, h - 1))
        return out
    def preds(self, word, h, real_lo=None, real_hi=None):
        lo = 0 if real_lo is None else real_lo
        hi = len(word) - 1 if real_hi is None else real_hi
        out = []
        if h + 1 <= hi:   # predecessor with head at h+1 moved left into T; theta(T->pred) is a right move of T
            out.append((self.Ecut(word, h - self.q + 1), word, h + 1))
        if h - 1 >= lo:
            out.append((1 / self.Ecut(word, h - self.q), word, h - 1))
        if word[h] == 0 and h - 1 >= lo and h + 1 <= hi:
            x, y = word[h - 1], word[h + 1]
            short = {(1, 1): 1, (1, 0): 0, (0, 1): 0}.get((x, y))
            if short is not None:
                w2 = word[:h - 1] + [short] + word[h + 2:]
                out.append((1 / self.Gv(w2, h - 1, (x, 0, y)), w2, h - 1))
        trips = [(0, 1, 0)] if word[h] == 0 else [(0, 1, 1), (1, 1, 0)]
        for tr in trips:
            w2 = word[:h] + list(tr) + word[h + 1:]
            out.append((self.Gv(word, h, tr), w2, h + 1))
        return out
    def row(self, word, h, real_lo=None, real_hi=None):
        W0 = self.W(word, h); Z0 = self.Z(word, h)
        s = Fr(0)
        for th, w2, h2 in self.moves(word, h, real_lo, real_hi):
            s += th * self.W(w2, h2) / W0
        for th, w2, h2 in self.preds(word, h, real_lo, real_hi):
            s += th * self.Z(w2, h2) / Z0
        return s / 2

PAD = 13
def formal_word(T, row):
    L = 2 * T['r'] + 3
    bits = [(row >> (L - 1 - i)) & 1 for i in range(L)]
    return [0] * PAD + bits + [0] * PAD, PAD + T['r'] + 1

def wordstr(bits, h):
    return ''.join(('[' + ('+' if b else '-') + ']') if i == h else ('+' if b else '-') for i, b in enumerate(bits))

def parse_opts(args):
    return dict(a.split('=') for a in args if '=' in a)

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'prep':
        opts = parse_opts(sys.argv[4:])
        T = load(sys.argv[2], opts)
        symmetry_report(T)
        r, k, q, dd = T['r'], T['k'], T['q'], T['d']
        with open(sys.argv[3], 'wb') as f:
            f.write(struct.pack('<4i', r, k, q, dd))
            f.write(struct.pack('<%dd' % len(T['U']), *[float(u) for u in T['U']]))
            f.write(struct.pack('<%dd' % len(T['F']), *[float(x) for x in T['F']]))
            f.write(struct.pack('<%dd' % len(T['E']), *[float(x) for x in T['E']]))
            nG = (1 << (2 * dd + 1)) * 8
            garr = [float(T['G'].get(i, Fr(0))) for i in range(nG)]
            f.write(struct.pack('<%dd' % nG, *garr))
        print('tables written', opts)
    elif cmd == 'exact':
        opts = parse_opts(sys.argv[5:])
        T = load(sys.argv[2], opts); Mc = Machine(T)
        Rref = Fr(sys.argv[4])
        rows = [int(l.split()[0]) for l in open(sys.argv[3]) if l.strip()]
        best = None; nbad = 0
        for row in rows:
            w, h = formal_word(T, row)
            v = Mc.row(w, h)
            if v > Rref:
                nbad += 1
            if best is None or v > best[0]:
                best = (v, row)
        print(f'exact re-evaluation of {len(rows)} flagged rows: {nbad} above Rref')
        if best:
            v, row = best
            w, h = formal_word(T, row)
            print(f'largest exact row value {float(v):.15f} at row {row} {wordstr(w[PAD:-PAD], h - PAD)}; equals Rref: {v == Rref}')
            print(f'R^2 = {float(v * v):.20f}; R^2 <= 12.2970: {v * v <= Fr(122970, 10000)}; R^2 <= 12.296933854: {v * v <= Fr(12296933854, 10**9)}')
            print('exact R =', v)
    elif cmd == 'sample':
        opts = parse_opts(sys.argv[4:])
        T = load(sys.argv[2], opts); Mc = Machine(T)
        rnd = random.Random(1390)
        for _ in range(int(sys.argv[3])):
            row = rnd.randrange(1 << (2 * T['r'] + 3))
            w, h = formal_word(T, row)
            print(row, repr(float(Mc.row(w, h))))
    elif cmd == 'tapes':
        opts = parse_opts(sys.argv[5:])
        T = load(sys.argv[2], opts); Mc = Machine(T)
        N = int(sys.argv[3]); maxlen = int(sys.argv[4]); rnd = random.Random(1391)
        r = T['r']; worst = Fr(0); nviol = 0; tested = 0
        while tested < N:
            # random walk from the root +[-]+ within tapes of at most maxlen letters
            word, h = [1, 0, 1], 1
            steps = rnd.randint(1, 4 * maxlen)
            for _ in range(steps):
                P = r + 2 + 2 * Mc.q
                pw = [0] * P + word + [0] * P
                mv = Mc.moves(pw, h + P, P, P + len(word) - 1)
                mv = [(a, b, c) for a, b, c in mv if len(b) - 2 * P <= maxlen]
                _, b, c = rnd.choice(mv)
                word, h = b[P:len(b) - P], c - P
            P = r + 2 + 2 * Mc.q
            pw = [0] * P + word + [0] * P
            actual = Mc.row(pw, h + P, P, P + len(word) - 1)   # full-tape weights, only real moves
            fw = pw[h + P - r - 1:h + P + r + 2]
            formal = Mc.row([0] * PAD + fw + [0] * PAD, PAD + r + 1)
            tested += 1
            if actual > formal:
                nviol += 1
            worst = max(worst, actual)
        print(f'actual tapes tested: {N} (length <= {maxlen}); actual row > formal row: {nviol}; largest actual row {float(worst):.12f}')
