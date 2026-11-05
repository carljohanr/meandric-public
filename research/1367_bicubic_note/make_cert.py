"""make_cert.py (1367): produce a flow certificate BICFLOW1 for the lower bound (float search; the
certificate itself is checked exactly by check_flow_a.py and check_flow_b.py).

Usage: python make_cert.py r sn sd iters out.bin
Reference law: symmetric two-state Markov chain, switch probability s = sn/sd (memory one).
File format BICFLOW1 (little endian):
  8 bytes magic 'BICFLOW1', int32 r, uint64 sn, sd, nA, nB, then nA uint64 X_A, then nB uint64 X_B.
  Windows are integers w = a + 2^r b + 2^(2r) c, where bit k-1 of a (of b) is the k-th symbol from the
  top of stack A (of stack B) and c is the parity.  X_A lists, in increasing w, the windows with
  (a mod 2) != c (pop A allowed); X_B likewise with (b mod 2) != c.  All entries are >= 1.
"""
import sys, os, struct, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fixpoint import run, entropy

r, sn, sd, iters, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
s = sn / sd
XA, XB, S, hist = run(r, 1, [s, s], iters, report=max(1, iters // 5))
selA, selB = S['popA_ok'], S['popB_ok']
xa, xb = XA[selA], XB[selB]
mx = max(xa.max(), xb.max())
ia = np.maximum(1, np.rint(xa / mx * 2.0 ** 60)).astype(np.uint64)
ib = np.maximum(1, np.rint(xb / mx * 2.0 ** 60)).astype(np.uint64)
print('smallest X / largest X: %.3e' % (min(xa.min(), xb.min()) / mx))
# float entropy of the rounded certificate
XA2 = np.zeros(S['N']); XB2 = np.zeros(S['N'])
XA2[selA] = ia.astype(float); XB2[selB] = ib.astype(float)
h, _, _ = entropy(S, XA2, XB2)
print('rounded certificate (float): h = %.12f, e^(2h) = %.10f' % (h, np.exp(2 * h)))
with open(out, 'wb') as f:
    f.write(b'BICFLOW1'); f.write(struct.pack('<i', r))
    f.write(struct.pack('<4Q', sn, sd, len(ia), len(ib)))
    f.write(ia.tobytes()); f.write(ib.tobytes())
print('wrote', out, os.path.getsize(out), 'bytes')
