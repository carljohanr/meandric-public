"""Keep accumulated free edits against a fixed exact parent, avoiding repeated normalisation denominators."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--previous',type=Path,required=True);ap.add_argument('--previous-basis',type=Path,required=True);ap.add_argument('--increment',type=Path,required=True);ap.add_argument('--basis',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.exists()
with np.load(a.previous) as z:old=z['effective_edit_coefficients'].copy()
with np.load(a.increment) as z:new=z['effective_edit_coefficients'].copy()
with np.load(a.previous_basis.with_name(a.previous_basis.name+'_metadata.npz')) as z:old_free=z['free'].copy();old_scales=list(zip(map(str,z['scale_numerators']),map(str,z['scale_denominators'])))
with np.load(a.basis.with_name(a.basis.name+'_metadata.npz')) as z:free=z['free'].copy();scales=list(zip(map(str,z['scale_numerators']),map(str,z['scale_denominators'])))
positions={int(c):i for i,c in enumerate(free)}
for c,value,scale in zip(old_free,old,old_scales):
 index=positions[int(c)];assert scales[index]==scale,('shared column normalization changed',int(c));new[index]+=value
meta=dict(status='ACCUMULATED_FLOAT_FREE_EDITS_FOR_FIXED_EXACT_PARENT',previous_sha256=hashlib.file_digest(a.previous.open('rb'),'sha256').hexdigest(),increment_sha256=hashlib.file_digest(a.increment.open('rb'),'sha256').hexdigest(),scope='Rounded stored sum; exact reconstruction and all original rows required before certification.')
np.savez(a.out,effective_edit_coefficients=new,metadata=np.array(json.dumps(meta)));print(json.dumps(dict(status=meta['status'],coordinates=len(new))))
