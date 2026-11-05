"""Exact rational coarse-parent reconstruction from free edits and head cuts.
Creates a diagnostic, not an entropy certificate. Signed intermediates require
--allow-signed and an independently verified positive completed fine flow.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import sys,time,json,pickle,argparse
sys.set_int_max_str_digits(0)
from pathlib import Path
import numpy as np
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'918_exact_backend/vendor'));sys.path.insert(0,str(ROOT/'1204_coarse_coordinates/vendor'))
from flint import fmpq,fmpz
from local_ratio_builder import build_local
from accurate_kernel import AccurateSymmetricHeadKernel
from stream_builder import windows
ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--allow-signed',action='store_true');ap.add_argument('--exact-base-parent',type=Path,required=True);ap.add_argument('--increment-witness',type=Path,required=True);ap.add_argument('--increment-basis-prefix',type=Path,required=True);ap.add_argument('--matrix-prefix',type=Path,required=True);ap.add_argument('--basis-prefix',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--export-physical',action='store_true');ap.add_argument('--recipe-only',action='store_true');a=ap.parse_args();assert not a.out.with_suffix('.json').exists();start=time.monotonic()
def emit(stage,**kw):
 rec=dict(stage=stage,seconds=time.monotonic()-start);rec.update(kw);print(json.dumps(rec),flush=True)
def fq(v):
 if isinstance(v,float):return fmpq(*v.as_integer_ratio())
 return fmpq(int(v.numerator),int(v.denominator))
with np.load(a.source,allow_pickle=False) as z:
 metadata=json.loads(str(z['metadata']));K=int(metadata['K']);r=int(metadata['r']);up={int(i):Q(v) for i,v in metadata['up'].items()};kappa=Q(metadata['kappa']);params=tuple(Q(v) for v in metadata.get('memory_parameters',[metadata.get('base_eta','1/8'),metadata['memory_eta']]));physical=z['physical'].copy();keys=z['keys'].copy()
s=build_local(K,r,up,kappa,params,exact_control=True);assert np.array_equal(keys,s.keys);j=physical/s.bw[s.vw]
recipe=pickle.load(a.matrix_prefix.with_name(a.matrix_prefix.name+'_matrix.pkl').open('rb'));recipe.pop('integer_rows',None)
if a.recipe_only:
 from types import SimpleNamespace
 from head_compatibility import compatibility
 C,heads,edits,reconstruct,td,ta=compatibility(s)
 with np.load(a.increment_basis_prefix.with_name(a.increment_basis_prefix.name+'_metadata.npz')) as bm:free=bm['free'].copy()
 ech=pickle.loads(a.matrix_prefix.with_name(a.matrix_prefix.name+'_echelon.pkl').read_bytes());pivots=ech['pivots'];priority=ech['priority'];del ech
 component=np.full(len(heads),-1,dtype=int);component[ta.roots]=np.arange(len(ta.roots))
 for level in ta.levels:component[level]=component[ta.parent[level]]
 head_index=np.full(len(s.keys),-1,dtype=int);head_index[heads]=np.arange(len(heads));partner=component[head_index[s.flowperm[heads[ta.roots]]]];representatives=np.minimum(np.arange(len(partner)),partner);unique,cut_orbit=np.unique(representatives,return_inverse=True)
 k=SimpleNamespace(free=free,pivots=pivots,back_order=sorted(pivots,key=lambda c:priority[c],reverse=True),B=SimpleNamespace(shape=(len(recipe['edits']),len(free))),heads=heads,edits=edits,edit_orbit=recipe['edit_orbit'],td=td,ta=ta,component=component,cut_orbit=cut_orbit,cut_count=len(unique))
else:k=AccurateSymmetricHeadKernel(s,a.matrix_prefix.with_name(a.matrix_prefix.name+'_matrix.pkl'),a.basis_prefix.with_suffix('.npz'),a.matrix_prefix.with_name(a.matrix_prefix.name+'_echelon.pkl'))
base_parent=pickle.load(a.exact_base_parent.open('rb'));rb=int(base_parent['r']);nb=2*rb+1;extra=r-rb;assert extra>=0
assert base_parent['K']==K and tuple(base_parent['parameters'])==params and base_parent['up']==up and base_parent['kappa']==kappa
base_lookup=np.full(K*(1<<nb)*8,-1,dtype=np.int32);base_lookup[base_parent['keys']]=np.arange(len(base_parent['keys']),dtype=np.int32);base_D=int(base_parent['denominator']);base_beta={};base_edits=[]
for fullcol in recipe['edits']:
 key=int(s.keys[int(fullcol)]);wk=key>>3;lev=wk>>(2*r+1);word=wk&((1<<(2*r+1))-1)
 for position in range(extra):lev+=1 if (word>>(2*r-position))&1 else -1
 wk0=(lev<<nb)|((word>>extra)&((1<<nb)-1));index=int(base_lookup[wk0*8+(key&7)]);assert index>=0
 if wk0 not in base_beta:
  bits=wk0&((1<<nb)-1);steps=tuple(1 if (bits>>q)&1 else -1 for q in range(nb-1,-1,-1));base_beta[wk0]=fq(s.law.exact.cylinder((lev,steps)))
 base_edits.append(fmpq(int(base_parent['numerators'][index]),base_D)/base_beta[wk0])
free_values=[fmpq(0)]*len(k.free)

with np.load(a.increment_basis_prefix.with_name(a.increment_basis_prefix.name+'_metadata.npz')) as z:
 selected_free=z['free'].copy();scales=[fmpq(int(n),int(d)) for n,d in zip(z['scale_numerators'],z['scale_denominators'])]
with np.load(a.increment_witness) as z:effective=z['effective_edit_coefficients'].copy()
assert len(effective)==len(scales)==len(selected_free)
# A COMMON dyadic denominator avoids an enormous lcm of unrelated column maxima.
precision=max(128,max(int(v.numerator).bit_length()-int(v.denominator).bit_length()+1 for v in scales)+96);common=1<<precision;position={int(c):i for i,c in enumerate(k.free)};rounding=fmpq(0)
for c,u,scale in zip(selected_free,effective,scales):
 target=fq(float(u))/scale;numerator,rem=divmod(int(target.numerator)*common,int(target.denominator));numerator+=int(2*rem>int(target.denominator));dyadic=fmpq(numerator,common);rounding=max(rounding,abs(dyadic-target)*scale);free_values[position[int(c)]]+=dyadic
assert rounding<fmpq(1,1<<94)
emit('exact_base_plus_dyadic_selected_increment',selected_columns=len(selected_free),common_dyadic_bits=precision,maximum_normalized_column_rounding=float(rounding))
del base_parent,base_lookup,base_beta

def backsolve(values,add_base=True):
 # exact forward substitution of the selected free-coordinate values through the retained echelon (no rational row cache)
 n=k.B.shape[0];x=[fmpq(0)]*n
 for c,v in zip(k.free,values):x[int(c)]=v
 for c in k.back_order:
  row=k.pivots[c];terms=[(int(v),x[int(i)]) for i,v in row.items() if i!=c and x[int(i)]]
  if terms:x[c]=-sum((v*t for v,t in terms),fmpq(0))/int(row[c])
 return [b+d for b,d in zip(base_edits,x)] if add_base else x

edits=backsolve(free_values);emit('edit_reconstruction',variables=len(edits))
beta=[fq(s.law.exact.cylinder((key>>(2*r+1),w))) for key,w in windows(K,2*r+1)];weights=[fmpq(0)]*len(edits)
for full,col in enumerate(k.edits):weights[int(k.edit_orbit[full])]+=beta[int(s.vw[col])]*int(s.db[col])
size=sum((w*v for w,v in zip(weights,edits)),fmpq(0))
if a.recipe_only:
 # Eliminate pivot variables from the exact size functional, without constructing
 # any nullspace basis columns. This is the transpose triangular substitution.
 size_cache=a.increment_basis_prefix.with_name(a.increment_basis_prefix.name+'_exact_size.npz')
 if size_cache.exists():
  import hashlib
  with np.load(size_cache) as cached:
   cached_meta=json.loads(str(cached['metadata']));assert cached_meta['basis_sha256']==hashlib.file_digest(a.increment_basis_prefix.with_suffix('.npz').open('rb'),'sha256').hexdigest();q=[fmpq(int(n),int(d))*scale for n,d,scale in zip(cached['numerators'],cached['denominators'],scales)]
 else:raise SystemExit('exact size cache required (normalize_transpose_size_cache.py)')
 assert sum((w*v for w,v in zip(weights,base_edits)),fmpq(0))==0
 assert size==sum((qv*v for qv,v in zip(q,free_values)),fmpq(0))
 adjoint=[fmpq(0)]*k.B.shape[0]
 for c,v in zip(k.free,q):adjoint[int(c)]=v

 scale_by_free={int(c):scale for c,scale in zip(selected_free,scales)};selected=max((position[int(c)] for c in selected_free),key=lambda i:abs(float(adjoint[int(k.free[i])]/scale_by_free[int(k.free[i])])))
else:selected=int(np.argmax(abs(k.size_gradient)))
if size:
 unit=[fmpq(0)]*len(k.free);unit[selected]=fmpq(1);direction=backsolve(unit,add_base=False);dsize=sum((w*v for w,v in zip(weights,direction)),fmpq(0));assert dsize!=0;correction=size/dsize;edits=[e-correction*d for e,d in zip(edits,direction)]
else:correction=fmpq(0)
assert sum((w*v for w,v in zip(weights,edits)),fmpq(0))==0
emit('exact_size_repaired',selected_free=int(k.free[selected]),correction_float=float(correction))
full_edits=[edits[int(o)] for o in k.edit_orbit];edit_index=np.full(len(s.keys),-1,dtype=np.int32);edit_index[k.edits]=np.arange(len(k.edits));forcing=[]
for row in s.exactrows:
 forcing.append(-sum((fq(v)*full_edits[int(edit_index[col])] for col,v in row.items() if edit_index[col]>=0),fmpq(0)))
def integrate(tree,values):
 x=[fmpq(0)]*tree.nodes
 for level in tree.levels:
  for v in level:x[int(v)]=x[int(tree.parent[v])]+int(tree.sign[v])*values[int(tree.edge[v])]
 return x
y=integrate(k.td,forcing);head=integrate(k.ta,y);emit('exact_heads_integrated',heads=len(head))
# Fit each additive head-cut constant to the desired parent, then reflect.
groups=k.cut_orbit[k.component];sums=np.zeros(k.cut_count);counts=np.bincount(groups,minlength=k.cut_count)
for i,col in enumerate(k.heads):sums[int(groups[i])]+=j[col]-float(head[i])
constants=[fq(float(v)) for v in sums/counts];head=[v+constants[int(groups[i])] for i,v in enumerate(head)];full=[fmpq(0)]*len(s.keys)
for col,v in zip(k.edits,full_edits):full[int(col)]=v
for col,v in zip(k.heads,head):full[int(col)]=v
full=[(full[i]+full[int(s.flowperm[i])])/2 for i in range(len(full))]
# Repair any negative exact reconstruction by adding nonnegative exact cycles.
# Choose the covering cycle with the least physical mass per repaired density.
cycle_repairs=[];added_mass=fmpq(0);nonpositive=[i for i,v in enumerate(full) if v<=0]
if nonpositive and not a.allow_signed:
 from local_cycle_family import kernel_local
 s.parameters=params;entries,descriptors=kernel_local(s,exact=True);columns=[{} for _ in descriptors];incidence={i:[] for i in nonpositive}
 for row,col,value in entries:columns[col][row]=columns[col].get(row,fmpq(0))+fq(value)
 masses=[]
 for col,rows in enumerate(columns):
  masses.append(sum((beta[int(s.vw[i])]*v for i,v in rows.items()),fmpq(0)))
  for i,v in rows.items():
   if i in incidence:incidence[i].append((col,v))
 ordered=sorted(nonpositive,key=lambda i:float(beta[int(s.vw[i])]*(-full[i])),reverse=True)
 for i in ordered:
  if full[i]>0:continue
  assert incidence[i],('uncovered',i)
  col,value=min(incidence[i],key=lambda cv:float(masses[cv[0]]/cv[1]));target=fq(max(1e-30,float(j[i])*1e-12));amount=(target-full[i])/value
  for row,coefficient in columns[col].items():full[row]+=amount*coefficient
  added_mass+=amount*masses[col];cycle_repairs.append((int(col),str(amount.numerator),str(amount.denominator)))
 emit('positive_cycle_repair',initial_nonpositive=len(nonpositive),cycles_added=len(cycle_repairs),physical_mass_added=float(added_mass))
 del entries,columns,incidence,masses

# Explicitly remove the remaining signed head current before final replay.
pre_drifts=[sum((beta[int(s.vw[i])]*int(drift[i])*v for i,v in enumerate(full)),fmpq(0)) for drift in (s.da,s.db)]
emit('pre_current_gauge',drift_a=float(pre_drifts[0]),drift_b=float(pre_drifts[1]),edit_sum=float(pre_drifts[0]+pre_drifts[1]))
assert pre_drifts[0]+pre_drifts[1]==0
current=[(-1 if int(key)&7==0 else 1 if int(key)&7==1 else 0) for key in s.keys]
current_drift=sum((beta[int(s.vw[i])]*int(s.da[i])*v for i,v in enumerate(current)),fmpq(0));current_correction=pre_drifts[0]/current_drift
assert abs(float(current_correction))<1e-8,('unexpected large current',float(current_correction))
full=[v-current_correction*c for v,c in zip(full,current)]
emit('exact_current_repaired',correction_float=float(current_correction))

# Exact original-row and both physical-drift replay. Do not rely only on recipe.
import math
common_density_denominator=fmpz(1)
for value in full:common_density_denominator=common_density_denominator.lcm(value.denominator)
density_numerators=[value.numerator*(common_density_denominator//value.denominator) for value in full]
minimum=min(full);neg=sum(v<=0 for v in full)
# Do not retain both the large rational reconstruction and integer replay vectors.
del full,full_edits,edits,y,head,forcing,base_edits
import gc
gc.collect()
emit('fraction_free_original_rows_start',common_density_bits=int(common_density_denominator).bit_length())
for rowno,row in enumerate(s.exactrows):
 rowden=1
 for value in row.values():rowden=math.lcm(rowden,int(value.denominator))
 residual=sum(int(value.numerator)*(rowden//int(value.denominator))*density_numerators[col] for col,value in row.items());assert residual==0,(rowno,'integer residual');s.exactrows[rowno]=None
 if rowno%10000==0:emit('exact_original_rows',rows=rowno)
# Multiply the two integer vectors once; avoid hundreds of thousands of rational
# normalisations and their large gcds. Original rows above remain density rows.
beta_denominator=fmpz(1)
for value in beta:beta_denominator=beta_denominator.lcm(value.denominator)
beta_numerators=[value.numerator*(beta_denominator//value.denominator) for value in beta]
physical_numerators=[value*beta_numerators[int(s.vw[i])] for i,value in enumerate(density_numerators)]
physical_denominator=common_density_denominator*beta_denominator
physical_mass=sum(physical_numerators,fmpz(0));assert physical_mass>0
for drift in (s.da,s.db):assert sum(int(d)*v for d,v in zip(drift,physical_numerators))==0
mass=fmpq(physical_mass,physical_denominator)
floats=np.asarray([int(v*mass.denominator)/int(common_density_denominator*mass.numerator) for v in density_numerators])
density_delta=floats-j/float(s.bw[s.vw]@j);physical_delta=float(np.sum(s.bw[s.vw]*abs(density_delta)));emit('exact_increment_proximity',maximum_density_difference=float(np.max(abs(density_delta))),physical_l1_difference=physical_delta);assert physical_delta<1e-8
del density_numerators,beta_numerators
gc.collect()
f=s.bw[s.vw]*floats;nu=np.bincount(s.vw,weights=f,minlength=len(s.bw));mu=float(np.exp(-2*f@np.log(f/nu[s.vw]))) if neg==0 else None
if a.export_physical:
 divisor=physical_mass
 for value in physical_numerators:
  divisor=divisor.gcd(value)
  if divisor==1:break
 D=int(physical_mass//divisor);numerators=[int(v//divisor) for v in physical_numerators]
 assert sum(numerators)==D and (a.allow_signed or min(numerators)>0)
 for drift in (s.da,s.db):assert sum(int(d)*v for d,v in zip(drift,numerators))==0
 physical_path=a.out.with_name(a.out.name+'_physical_certificate.pkl')
 import hashlib
 digest=hashlib.sha256()
 class CheckedWriter:
  def __init__(self,target):self.target=target
  def write(self,data):
   written=self.target.write(data);assert written==memoryview(data).nbytes;digest.update(data);return written
 with physical_path.open('xb') as fcert:
  pickle.dump(dict(K=K,r=r,up=up,kappa=kappa,parameters=params,keys=s.keys,denominator=D,numerators=numerators),CheckedWriter(fcert),protocol=5);fcert.flush();os.fsync(fcert.fileno())
 intended_digest=digest.hexdigest();observed_digest=hashlib.file_digest(physical_path.open('rb'),'sha256').hexdigest();assert observed_digest==intended_digest,('export read-back mismatch',intended_digest,observed_digest)
 emit('physical_integer_parent_exported',common_denominator_bits=D.bit_length(),bytes=physical_path.stat().st_size,sha256=observed_digest,byte_stream_hash_verified=True)
result=dict(signed_intermediate=a.allow_signed,exact_base_parent=str(a.exact_base_parent),increment_witness=str(a.increment_witness),increment_normalized_rounding=float(rounding),physical_l1_difference=physical_delta,recipe_only=a.recipe_only,status='EXACT_RATIONAL_PARENT_DIAGNOSTIC',K=K,r=r,up=metadata['up'],kappa=str(kappa),memory_order=len(params)+1,memory_parameters=list(map(str,params)),memory_eta=str(params[-1]),base_eta=str(params[0]),all_original_rows_exact=True,both_drifts_exact=True,nonpositive_flows=neg,positive_cycle_repairs=len(cycle_repairs),cycle_repair_added_mass=float(added_mass),minimum_density_float=float(minimum),mu_float_not_certified=mu,source=str(a.source),seconds=time.monotonic()-start)
a.out.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');np.savez(a.out.with_suffix('.npz'),keys=s.keys,physical=f,metadata=np.array(json.dumps(result)))
# Compact exact recipe: initial free edits, one exact correction and cut constants.
compact=dict(exact_base_parent=str(a.exact_base_parent),increment_basis_prefix=str(a.increment_basis_prefix),cycle_repairs=cycle_repairs,source=str(a.source),free_numerators=[str(v.numerator) for v in free_values],free_denominators=[str(v.denominator) for v in free_values],correction_free_index=selected,correction_numerator=str(correction.numerator),correction_denominator=str(correction.denominator),cut_numerators=[str(v.numerator) for v in constants],cut_denominators=[str(v.denominator) for v in constants],normalization_numerator=str(mass.numerator),normalization_denominator=str(mass.denominator),metadata=result)
a.out.with_name(a.out.name+'_exact_recipe.json').write_text(json.dumps(compact)+'\n');emit('finished',**result)
