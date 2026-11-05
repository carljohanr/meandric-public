"""Write a verify.py-compatible candidate JSON from a Model3 (no refinement) and x."""
import json,hashlib
import numpy as np
from fractions import Fraction as Q
def save_candidate(m,x,path,extra=None):
    assert not m.selU and not m.selE and not m.selG and not m.selF
    Uo=np.maximum(1,np.rint(np.exp(x[:m.nu]-max(x[:m.nu]))*2**48)).astype(np.int64)
    def quant(v):return str(Q(max(1,round(float(np.exp(v))*2**36)),2**36))
    cand=dict(status='NUMERICAL_CANDIDATE_REQUIRES_EXACT_REPLAY',r=m.rU,k=m.k,q=m.q,d=m.d,
        U=Uo[m.lookup].tolist(),F=[quant(v) for v in x[m.oF:m.oF+(1<<m.k)]],E=[quant(v) for v in x[m.oE:m.oE+m.ne]],G=[quant(v) for v in x[m.oG:m.oG+m.nd]],
        field_lookup=m.field,edit_lookup={str(a):b for a,b in m.emap.items()})
    if extra:cand.update(extra)
    with open(path,'w') as f:json.dump(cand,f);f.write('\n')
