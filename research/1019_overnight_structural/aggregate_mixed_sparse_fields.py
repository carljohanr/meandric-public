"""Collapse same-radius short/six/eight stages into one literal field witness.
Every rounded coefficient is reinterpreted exactly by the independent verifier.
"""
from pathlib import Path
import argparse,json,hashlib,numpy as np,time
ap=argparse.ArgumentParser();ap.add_argument('--short',type=Path,required=True);ap.add_argument('--six',type=Path,required=True);ap.add_argument('--eight',type=Path,required=True);ap.add_argument('--joint',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.exists();start=time.monotonic();paths=[a.short,a.six,a.eight,a.joint];data=[];metadata=[]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
for p in paths:
 with np.load(p) as z:data.append({k:z[k].copy() for k in ('descriptors','coefficients','changed_keys','density_increment')});metadata.append(json.loads(str(z['metadata'])))
for prev,meta in zip(paths,metadata[1:]):assert Path(meta['source']).resolve()==prev.resolve() and meta['source_sha256']==sha(prev)
assert len({m['r'] for m in metadata})==len({m['K'] for m in metadata})==1 and metadata[2]['family']=='eight_step_excursion_literal_cycles' and metadata[3]['family']=='mixed_literal_cycles'
coef={};keys=data[-1]['changed_keys'];delta=np.zeros(len(keys))
for index,z in enumerate(data):
 offset=4 if index==2 else 0
 for (kind,key),value in zip(z['descriptors'],z['coefficients']):d=(int(kind)+offset,int(key));coef[d]=coef.get(d,0.)+float(value)
 ii=np.searchsorted(keys,z['changed_keys']);assert np.all(ii<len(keys)) and np.array_equal(keys[ii],z['changed_keys']);delta[ii]+=z['density_increment']
desc=sorted(coef);meta=dict(metadata[-1],source=metadata[0]['source'],source_sha256=metadata[0]['source_sha256'],initial_mu=metadata[0]['initial_mu'],directions=len(desc),changed_flows=len(keys),changed_windows=len(np.unique(keys>>3)),minimum_changed_density=min(m['minimum_changed_density'] for m in metadata),aggregated_sources=[dict(path=str(p),sha256=sha(p)) for p in paths],aggregation_note='Binary-float sums require independent exact path replay; preceding certificates are not assumed to certify these rounded sums.',seconds=time.monotonic()-start)
np.savez(a.out,descriptors=np.asarray(desc,dtype=np.int64),coefficients=np.asarray([coef[d] for d in desc]),changed_keys=keys,density_increment=delta,metadata=np.array(json.dumps(meta)));a.out.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta))
