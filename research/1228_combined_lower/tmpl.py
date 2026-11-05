"""Template generators for the 1211 trial: maximal-context cycles and source-fitting diamonds.

A template is a closed signed walk on literal tape configurations (word, head).  Its column is
Gamma(f) = sum_e f(e) gamma(e), gamma(e) = (1/beta(source word)) * indicator of every radius-r window completing
the source configuration, with the edge label.  Only SOURCES must fit the window (Lemma 1 of 1205).
Cycle templates: all states are sources.  Diamond templates: two sources a, c; targets b, d need not fit.

Coefficient of the edge column with source word w_t is  sign * beta(anchor)/beta(w_t)  (exact rational,
converted to float in the screen/optimizer; the verifier recomputes it exactly).
"""
import sys
from pathlib import Path
PROD = Path(__file__).resolve().parent.parent / '1019_overnight_structural'
sys.path.insert(0, str(PROD))
from sparse_fine_cycles import *   # noqa: F401,F403  (literal, windows, code, LABELS, itertools, Q, np, ...)

BASE_KIND = 40          # template index = kind - BASE_KIND
D_LEN = {'L': 0, 'R': 0, 'S1': 2, 'S2': 2, 'S3': 2, 'M': -2}
D_HEAD = {'L': -1, 'R': 1, 'S1': 1, 'S2': 1, 'S3': 1, 'M': -1}

# four-step cycles (cyclic label sequences, S first) and the 14 eight-step classes of 1205 with positive rank gain
FOUR = (('S1', 'L', 'M', 'R'), ('S2', 'L', 'M', 'R'), ('S1', 'R', 'M', 'L'), ('S3', 'R', 'M', 'L'))
EIGHT_1205 = (
    'L L M R R M S1 S3', 'L L M R R M S2 S3', 'L L M S1 S2 R R M', 'L L M S3 S2 R R M', 'L L M M R R S1 S2',
    'L L S1 S3 R R M M', 'L L M M S1 R R S2', 'L L M M S3 R R S1', 'L L S1 R R M M S2', 'L L S3 R R M M S1',
    'L L M M S1 R R S1', 'L L M M S3 R R S2', 'L L S1 R R M M S1', 'L L S3 R R M M S2')
MIRROR_LAB = {'L': 'R', 'R': 'L', 'S1': 'S1', 'S2': 'S3', 'S3': 'S2', 'M': 'M'}


def canon(seq):
    seq = tuple(seq)
    return min(seq[i:] + seq[:i] for i in range(len(seq)))


def trajectory(labels):
    """(dlen, dhead) of the state BEFORE each label, relative to state 0."""
    out = [(0, 0)]
    for lab in labels[:-1]:
        out.append((out[-1][0] + D_LEN[lab], out[-1][1] + D_HEAD[lab]))
    return out


def max_context(labels, r):
    """Largest (m, h0) with all states fitting radius r.  m = length of state 0."""
    tr = trajectory(labels)
    mh = max(dh for _, dh in tr)
    mrt = max(dl - dh for dl, dh in tr)
    h0 = r - mh
    m = h0 + 1 + r - mrt
    return m, h0


class Template:
    def __init__(self, name, kind, labels, m, h0, dtype='cycle', side=None, split=None):
        self.name = name; self.kind = kind; self.labels = tuple(labels); self.m = m; self.h0 = h0
        self.dtype = dtype; self.side = side; self.split = split

    def as_dict(self):
        return dict(name=self.name, kind=self.kind, labels=list(self.labels), m=self.m, h0=self.h0, dtype=self.dtype,
                    side=self.side, split=self.split)


def cycle_template(name, kind, labels, r, m=None, h0=None):
    mm, hh = max_context(labels, r)
    return Template(name, kind, labels, mm if m is None else m, hh if h0 is None else h0)


PROD6 = (('L', 'S1', 'S3', 'R', 'M', 'M'), ('L', 'S3', 'S3', 'R', 'M', 'M'), ('R', 'S1', 'S2', 'L', 'M', 'M'), ('R', 'S2', 'S2', 'L', 'M', 'M'))
PROD8 = (('L', 'L', 'S1', 'R', 'R', 'R', 'M', 'L'), ('L', 'L', 'S3', 'R', 'R', 'R', 'M', 'L'), ('L', 'S1', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S1', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S1', 'S3', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S3', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'S3', 'R', 'R', 'M', 'M', 'L'), ('R', 'R', 'S1', 'L', 'L', 'L', 'M', 'R'), ('R', 'R', 'S2', 'L', 'L', 'L', 'M', 'R'), ('R', 'S1', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S1', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S1', 'S3', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S2', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'S3', 'L', 'L', 'M', 'M', 'R'))


def build_registry(r, groups=('c4', 'e8', 'dia'), mirror_close=True):
    """Templates of radius r.  Kinds are FIXED (order p6, p8, c4, e8, dia) so descriptors are stable across runs.
    Groups: c4 e8 dia (new), p6 p8 (production patterns re-run at maximal context)."""
    full = []
    for lab in PROD6:
        full.append(('p6', 'p6:' + ' '.join(lab), lab))
    seen = []
    for lab in PROD8:
        c = canon(lab)
        if c not in seen:
            seen.append(c)
            full.append(('p8', 'p8:' + ' '.join(c), c))
    for lab in FOUR:
        full.append(('c4', 'c4:' + ' '.join(lab), lab))
    classes = []
    for s in EIGHT_1205:
        c = canon(s.split())
        if c not in classes:
            classes.append(c)
        if mirror_close:
            mc = canon([MIRROR_LAB[x] for x in c])
            if mc not in classes:
                classes.append(mc)
    for c in classes:
        full.append(('e8', 'e8:' + ' '.join(c), c))
    out = []
    for idx, (g, name, lab) in enumerate(full):
        if g in groups:
            out.append(cycle_template(name, BASE_KIND + idx, lab, r))
    if 'dia' in groups:
        out.extend(diamond_templates(r, BASE_KIND + len(full)))
    if 'p4' in groups:   # production four-step contexts (m = 2r-3, centered head), for matched-budget controls
        for i, lab in enumerate(FOUR):
            out.append(cycle_template('p4:' + ' '.join(lab), BASE_KIND + len(full) + 4 + i, lab, r, m=2 * r - 3, h0=r - 2))
    if 'q6' in groups:   # production six-step contexts exactly (m = 2r-5, centered head, rotations as in 1019)
        for i, lab in enumerate(PROD6):
            out.append(cycle_template('q6:' + ' '.join(lab), BASE_KIND + len(full) + 8 + i, lab, r, m=2 * r - 5, h0=(2 * r - 5) // 2))
    if 'q8' in groups:   # production eight-step contexts exactly
        for i, lab in enumerate(PROD8):
            out.append(cycle_template('q8:' + ' '.join(lab), BASE_KIND + len(full) + 12 + i, lab, r, m=2 * r - 5, h0=(2 * r - 5) // 2))
    return out


def diamond_templates(r, kind0):
    out = []
    for side, split in (('L', 'S1'), ('L', 'S2'), ('R', 'S1'), ('R', 'S3')):
        h0 = r if side == 'L' else r - 2
        out.append(Template(f'dia:{side}:{split}', kind0 + len(out), (side, split), 2 * r - 1, h0,
                            dtype='diamond', side=side, split=split))
    return out


def diff_window(wa, wt, collar=3):
    """Local sub-words of wa, wt (aligned at common prefix/suffix) containing all differences plus collars."""
    p = 0
    n = min(len(wa), len(wt))
    while p < n and wa[p] == wt[p]:
        p += 1
    s = 0
    while s < n - p and wa[len(wa) - 1 - s] == wt[len(wt) - 1 - s]:
        s += 1
    lo = max(p - collar, 0)
    return lo, wa[lo:len(wa) - s + collar], wt[lo:len(wt) - s + collar]


class TemplateGenerator:
    """Descriptor (kind, key): key = (lev << m) | code(anchor word of length m)."""
    def __init__(self, K, r, law, parameters, groups=('c4', 'e8', 'dia'), templates=None):
        self.K = K; self.r = r; self.n = 2 * r + 1; self.law = law; self.parameters = parameters
        self.memory = len(parameters) + 1
        self.templates = templates if templates is not None else build_registry(r, groups)
        self.by_kind = {t.kind: t for t in self.templates}
        self.ratios = {}

    # ---- nodes and edges of one instance ----
    def instance(self, t, lev, anchor):
        """Return list of (word, head, label, sign, source_index_for_ratio) or [] if the instance is illegal."""
        K = self.K
        a = (anchor, t.h0)
        if t.dtype == 'cycle':
            nodes = []; w, h = a
            for lab in t.labels:
                nodes.append((w, h, lab, 1))
                try:
                    w, h = literal.literal_move(w, h, lab)
                except ValueError:
                    return []
                if not literal.admissible(lev, w, K) or w[:3] != anchor[:3] or w[-3:] != anchor[-3:]:
                    return []
            if (w, h) != a:
                return []
            return nodes
        # diamond
        side, split = t.side, t.split
        back = 'R' if side == 'L' else 'L'
        try:
            b = literal.literal_move(anchor, t.h0, side)
            d = literal.literal_move(anchor, t.h0, split)
        except ValueError:
            return []
        c = (d[0], d[1] - 1 if side == 'L' else d[1] + 1)
        try:
            mc = literal.literal_move(c[0], c[1], 'M')
        except (ValueError, AssertionError):
            return []
        if mc != b:
            return []
        for (w, h) in (b, c, d):
            if not literal.admissible(lev, w, K) or w[:3] != anchor[:3] or w[-3:] != anchor[-3:]:
                return []
        return [(anchor, t.h0, side, 1), (anchor, t.h0, split, -1), (c[0], c[1], 'M', -1), (c[0], c[1], back, 1)]

    def ratio(self, lev, wa, wt, exact=False):
        if wa == wt:
            return Q(1)
        lo, la, lt = diff_window(wa, wt)
        start = lev + sum(wa[:lo])
        key = (start, la, lt)
        if key not in self.ratios:
            auth = getattr(self.law, 'exact', self.law)
            self.ratios[key] = auth.cylinder((start, la)) / auth.cylinder((start, lt))
        return self.ratios[key]

    def entries(self, descriptor, exact=False, word=None):
        kind, key = map(int, descriptor)
        t = self.by_kind[kind]; m = t.m
        lev = key >> m
        anchor = word if word is not None else literal.decode(key & ((1 << m) - 1), m)
        if not literal.admissible(lev, anchor, self.K):
            return []
        nodes = self.instance(t, lev, anchor)
        if not nodes:
            return []
        r, n = self.r, self.n
        out = []
        for (w, h, lab, sign) in nodes:
            left = r - h; right = r - (len(w) - 1 - h)
            if left < 0 or right < 0:
                return []
            value = sign * self.ratio(lev, anchor, w)
            if not exact:
                value = float(value)
            for exterior in itertools.product((1, -1), repeat=left + right):
                prefix = exterior[:left]
                full = prefix + w + exterior[left:]
                start = lev - sum(prefix)
                if literal.admissible(start, full, self.K):
                    out.append(((((start << n) | code(full)) << 3) | LABELS.index(lab), value))
        return out
