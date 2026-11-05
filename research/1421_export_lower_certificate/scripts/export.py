#!/usr/bin/env python3
"""1421 export: write the Paper B lower-bound certificate in documented plain formats (see ../SPEC.md).

Reads the project's own files (the pickled radius-6 parent and the five npz stages) with the project's own decoding
rules (940 `literal`, 1019 `MemoryThreeLaw`, the descriptor rules of the 1019/1228 verifiers), and writes:
  export/law.json                 reference law: parameters and exact tables
  export/parent_r6.bin            radius-6 signed parent, binary (format IWPAR001, SPEC section 4)
  export/parent_r6_index.tsv      the same parent, one text line per flow, with float approximations only
  export/dirs_<s>_<name>.tsv      closed-walk directions of stage s, one line per direction (SPEC section 5)
  export/stages.json              stage list, counts, recorded certificate values, native file hashes
Run from anywhere:  python export.py   (uses ../export relative to this script).
"""
import sys, json, hashlib, pickle, time, itertools, struct
from pathlib import Path
from fractions import Fraction as Q
import numpy as np

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
RES = FOLDER.parent                      # meandric/research
OUT = FOLDER / 'export'
sys.path[:0] = [str(RES / '1019_overnight_structural'), str(RES / '910_band_bound_push/vendor')]
sys.set_int_max_str_digits(0)
import importlib.util
spec = importlib.util.spec_from_file_location('literal940', RES / '940_integer_literal_verifier/verify_sorted.py')
literal = importlib.util.module_from_spec(spec); spec.loader.exec_module(literal)
from memory_three import MemoryThreeLaw

A7 = RES / '1204_coarse_coordinates/data/alt7'
C28 = RES / '1228_combined_lower'
PARENT = A7 / 'exact_parent_physical_certificate.pkl'
STAGES = [
    # (index, short name, npz, certificate json, radius)
    (0, 'base_r8', A7 / 'combined_cycle_r8.npz', A7 / 'r8_certificate.json', 8),
    (1, 'sparse_r9', A7 / 'sparse_r9_candidate.npz', A7 / 'sparse_r9_certificate.json', 9),
    (2, 'aggregated_r10', A7 / 'aggregated_mixed_r10_candidate.npz', A7 / 'aggregated_mixed_r10_certificate.json', 10),
    (3, 'increment1_r10', C28 / 'runs/inc1_131072.npz', C28 / 'inc1_131072_certificate.json', 10),
    (4, 'increment2_r10', C28 / 'runs/inc2_131072.npz', C28 / 'inc2_131072_certificate.json', 10),
]
SIX = (('L','S1','S3','R','M','M'),('L','S3','S3','R','M','M'),('R','S1','S2','L','M','M'),('R','S2','S2','L','M','M'))
EIGHT = (('L','L','S1','R','R','R','M','L'),('L','L','S3','R','R','R','M','L'),('L','S1','R','R','R','M','L','L'),
         ('L','S1','S2','R','R','M','M','L'),('L','S1','S3','R','R','M','M','L'),('L','S3','R','R','R','M','L','L'),
         ('L','S3','S2','R','R','M','M','L'),('L','S3','S3','R','R','M','M','L'),('R','R','S1','L','L','L','M','R'),
         ('R','R','S2','L','L','L','M','R'),('R','S1','L','L','L','M','R','R'),('R','S1','S2','L','L','M','M','R'),
         ('R','S1','S3','L','L','M','M','R'),('R','S2','L','L','L','M','R','R'),('R','S2','S2','L','L','M','M','R'),
         ('R','S2','S3','L','L','M','M','R'))
LET = {1: '+', -1: '-'}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 23), b''): h.update(b)
    return h.hexdigest()


def word_str(w): return ''.join(LET[x] for x in w)


def frac(q): return f'{q.numerator}/{q.denominator}' if q.denominator != 1 else str(q.numerator)


def decode_direction(kind, key, radius, family, table):
    """Return (level, anchor word, head index, labels) exactly as the project's verifiers decode them."""
    n = 2 * radius + 1
    if family == 'template_trial_cycles':
        t = table[kind]; assert t['dtype'] == 'cycle'
        m, h = int(t['m']), int(t['h0'])
        return key >> m, literal.decode(key & ((1 << m) - 1), m), h, tuple(t['labels'])
    if kind == 0:                         # head cut: tape of 2r letters, head at index r, walk L then R
        lev = key >> (n - 1); w = literal.decode(key & ((1 << (n - 1)) - 1), n - 1)
        return lev, w, radius, ('L', 'R')
    if kind in (1, 2):                    # four-move edit cycle, anchor of n-4 letters, head at its middle
        m = n - 4; lev = key >> m; w = literal.decode(key & ((1 << m) - 1), m); h = m // 2
        split = 'S1' if w[h] == 1 else ('S2' if kind == 1 else 'S3')
        labels = (split, 'L', 'M', 'R') if kind == 1 else (split, 'R', 'M', 'L')
        return lev, w, h, labels
    assert family == 'mixed_literal_cycles' and 3 <= kind < 3 + len(SIX) + len(EIGHT)
    m = n - 6; lev = key >> m; w = literal.decode(key & ((1 << m) - 1), m)
    return lev, w, m // 2, (SIX + EIGHT)[kind - 3]


def walk_check(lev, anchor, head, labels, radius, K):
    """Replay the walk; assert closure, protection (level, 3+3 end letters, letter sum), interior head, fit, zero drift."""
    word, h = anchor, head; d = [0, 0]; maxlen = len(anchor)
    for lab in labels:
        assert literal.admissible(lev, word, K)
        assert word[:3] == anchor[:3] and word[-3:] == anchor[-3:] and sum(word) == sum(anchor)
        assert 0 < h < len(word) - 1
        assert radius - h >= 0 and radius - (len(word) - 1 - h) >= 0, 'source tape does not fit the stage radius'
        dd = literal.DRIFT[lab]; d[0] += dd[0]; d[1] += dd[1]
        word, h = literal.literal_move(word, h, lab); maxlen = max(maxlen, len(word))
    assert (word, h) == (anchor, head) and d == [0, 0]
    return maxlen


def export_law(K, up, kappa, eta, zeta):
    law = MemoryThreeLaw(K, up, kappa, eta, zeta)
    obj = dict(format='IW-LAW-1', K=K, memory=3, up=[frac(up[i]) for i in range(K)], kappa=frac(kappa),
               eta=frac(eta), zeta=frac(zeta),
               note='heights 0..K-1; tables keyed by consecutive heights; probabilities are exact fractions',
               rho=[[*k, frac(v)] for k, v in sorted(law.rho.items())],
               pi=[[*k, frac(v)] for k, v in sorted(law.pi.items())],
               quartet=[[*k, frac(v)] for k, v in sorted(law.quartet.items())],
               quintet=[[*k, frac(v)] for k, v in sorted(law.quintet.items())])
    (OUT / 'law.json').write_text(json.dumps(obj, indent=1) + '\n')
    return law


def export_parent(c, law):
    K, r = c['K'], c['r']; n = 2 * r + 1
    D = c['denominator']; N = c['numerators']; keys = c['keys']
    # canonical order: level ascending, word lexicographic with + before -, labels in order L R S1 S2 S3 M
    canon = []
    for lev, w in literal.windows(K, n):
        for lab in literal.legal_labels(lev, w, r, K):
            canon.append((((lev << n) | literal.code(w)) << 3) | literal.LABELS.index(lab))
    assert np.array_equal(np.asarray(canon, dtype=np.uint64), keys.astype(np.uint64)), 'pickle not in canonical order'
    assert sum(N) == D
    with open(OUT / 'parent_r6.bin', 'wb') as f, open(OUT / 'parent_r6_index.tsv', 'w') as t:
        db = D.to_bytes((D.bit_length() + 7) // 8, 'little')
        f.write(b'IWPAR001'); f.write(struct.pack('<IIII', K, r, len(N), len(db))); f.write(db)
        t.write('# 1421 radius-6 parent index (float approximations only; exact values are in parent_r6.bin)\n')
        t.write('# columns: index level word label sign numerator_bits flow_approx density_approx\n')
        for i, (key, x) in enumerate(zip(keys, N)):
            key = int(key); lab = key & 7; wk = key >> 3; lev = wk >> n; code = wk & ((1 << n) - 1)
            mag = abs(x); mb = mag.to_bytes((mag.bit_length() + 7) // 8, 'little')
            f.write(struct.pack('<BBBBII', lev, lab, 1 if x < 0 else 0, 0, code, len(mb))); f.write(mb)
            w = literal.decode(code, n); flow = float(Q(x, D)); dens = flow / float(law.cylinder((lev, w)))
            t.write(f'{i}\t{lev}\t{word_str(w)}\t{literal.LABELS[lab]}\t{"-" if x < 0 else "+"}\t{mag.bit_length()}'
                    f'\t{flow:.17g}\t{dens:.17g}\n')
    return dict(K=K, r=r, flows=len(N), denominator_bits=D.bit_length(), negative=sum(x < 0 for x in N))


def export_stage(idx, name, npz, radius, K, law):
    z = np.load(npz, allow_pickle=False)
    desc = z['descriptors']; coef = z['coefficients']; meta = json.loads(str(z['metadata']))
    assert int(meta['r']) == radius and int(meta['K']) == K
    family = meta.get('family', 'four_step_and_head_cut')
    table = {int(t['kind']): t for t in meta.get('templates', [])}
    fn = f'dirs_{idx}_{name}.tsv'; counts = {}; maxlen = 0
    with open(OUT / fn, 'w') as f:
        f.write(f'# 1421 closed-walk directions, stage {idx} ({name}), stage radius {radius}, K {K}\n')
        f.write('# columns: index level anchor head moves coef_hex coef_num coef_pow2 native_kind native_key\n')
        f.write('# coefficient c = coef_num / 2^coef_pow2 exactly (= float64 coef_hex); see SPEC.md section 5\n')
        for i, ((kind, key), v) in enumerate(zip(desc, coef)):
            kind, key, v = int(kind), int(key), float(v)
            assert np.isfinite(v)
            lev, w, h, labels = decode_direction(kind, key, radius, family, table)
            maxlen = max(maxlen, walk_check(lev, w, h, labels, radius, K))
            num, den = v.as_integer_ratio(); k = den.bit_length() - 1; assert den == 1 << k
            f.write(f'{i}\t{lev}\t{word_str(w)}\t{h}\t{",".join(labels)}\t{v.hex()}\t{num}\t{k}\t{kind}\t{key}\n')
            sig = (len(w), h, ','.join(labels)); counts[sig] = counts.get(sig, 0) + 1
    return dict(index=idx, name=name, radius=radius, directions=int(len(desc)), file=fn,
                walk_types=[dict(anchor_length=a, head=h, moves=mv, count=c) for (a, h, mv), c in sorted(counts.items())],
                max_tape_length=maxlen, native_family=family)


def main():
    t0 = time.monotonic(); OUT.mkdir(exist_ok=True)
    c = pickle.load(open(PARENT, 'rb'))
    K = c['K']; up = c['up']; kappa = c['kappa']; eta, zeta = c['parameters']
    law = export_law(K, up, kappa, eta, zeta)
    pinfo = export_parent(c, law); del c
    print(json.dumps(dict(stage='parent', seconds=time.monotonic() - t0, **pinfo)), flush=True)
    stages = []
    for idx, name, npz, cert, radius in STAGES:
        info = export_stage(idx, name, npz, radius, K, law)
        cj = json.loads(Path(cert).read_text())
        info['native_npz'] = str(npz.relative_to(RES)); info['native_npz_sha256'] = sha(npz)
        info['project_certificate'] = str(Path(cert).relative_to(RES)); info['project_certificate_sha256'] = sha(cert)
        info['recorded'] = {k: cj[k] for k in ('mu_lower_decimal', 'mu_lower_fraction', 'entropy_lower', 'mass_interval',
                            'changed_flows', 'changed_windows', 'flows_in_changed_windows_checked', 'windows', 'flows',
                            'minimum_density_lower', 'minimum_checked_density_lower') if k in cj}
        stages.append(info)
        print(json.dumps(dict(stage=name, seconds=time.monotonic() - t0, directions=info['directions'])), flush=True)
    obj = dict(format='IW-STAGES-1', K=K, final_radius=10, parent=dict(file='parent_r6.bin', index='parent_r6_index.tsv',
               native_pickle=str(PARENT.relative_to(RES)), native_pickle_sha256=sha(PARENT), **pinfo),
               law='law.json', stages=stages, total_directions=sum(s['directions'] for s in stages))
    (OUT / 'stages.json').write_text(json.dumps(obj, indent=1) + '\n')
    print(json.dumps(dict(stage='done', seconds=time.monotonic() - t0, total=obj['total_directions'])), flush=True)


if __name__ == '__main__':
    main()
