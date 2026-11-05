"""Packed geometry and literal cycles with exact short-ratio coefficients.
No fine stationarity matrix is constructed. Small exact row controls required.
"""
import sys,itertools,time,json
from pathlib import Path
from array import array
from types import SimpleNamespace
from fractions import Fraction as Q
import numpy as np
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parent.parent
sys.path[:0]=[str(ROOT/'1008_lower_122'),str(ROOT/'1017_transport_flows')]
from cycle_kernel import literal,windows,code,LABELS
from memory_order import FloatOrderLaw
from memory_two import FloatMemoryTwo
from float_system import FloatLaw

def geometry_witness(source,r):
 with np.load(source,allow_pickle=False) as z:
  metadata=json.loads(str(z['metadata']));K=int(metadata['K']);r0=int(metadata['r']);up={int(i):Q(v) for i,v in metadata['up'].items()};kappa=Q(metadata['kappa'])
  if 'memory_parameters' in metadata:parameters=tuple(map(Q,metadata['memory_parameters']))
  elif metadata.get('memory_order')==3:parameters=(Q(metadata.get('base_eta','1/8')),Q(metadata['memory_eta']))
  elif metadata.get('memory_eta') is not None:parameters=(Q(metadata['memory_eta']),)
  else:parameters=()
  law=FloatLaw(K,up,kappa) if not parameters else (FloatMemoryTwo(K,up,kappa,parameters[0]) if len(parameters)==1 else FloatOrderLaw(K,up,kappa,parameters))
  if 'memory_profile' in metadata:
   from mixed_profile_law import make_law
   law=make_law(K,up,kappa,parameters,{int(g):Q(v) for g,v in metadata['memory_profile'].items()})
  n0=2*r0+1;lookup=np.full(K*(1<<n0)*8,-1,dtype=np.int32);lookup[z['keys']]=np.arange(len(z['keys']),dtype=np.int32)
  old_weights=[]
  for key,w in windows(K,n0):
   weight=law.cylinder((key>>n0,w))
   old_weights.extend([weight]*len(literal.legal_labels(key>>n0,w,r0,K)))
  density=z['physical']/np.asarray(old_weights)
 n=2*r+1;extra=r-r0;assert extra>=0
 keys=array('I');vw=array('i');bw=array('d');base=array('d');da=array('b');db=array('b')
 for wi,(wk,w) in enumerate(windows(K,n)):
  lev=wk>>n;bw.append(law.cylinder((lev,w)));oldword=w[extra:n-extra] if extra else w;oldlev=lev+sum(w[:extra]);oldkey=((oldlev<<n0)|code(oldword))*8
  for label in literal.legal_labels(lev,w,r,K):
   li=LABELS.index(label);index=lookup[oldkey+li];assert index>=0
   keys.append(wk*8+li);vw.append(wi);base.append(density[index]);d=literal.DRIFT[label];da.append(d[0]);db.append(d[1])
 return SimpleNamespace(K=K,r=r,keys=np.asarray(keys),vw=np.asarray(vw),bw=np.asarray(bw),base=np.asarray(base),da=np.asarray(da),db=np.asarray(db),law=law,parameters=parameters,source_metadata=metadata)

def kernel_local(s,exact=False,include_heads=True):
 K,r=s.K,s.r;n=2*r+1;memory=len(s.parameters)+1;assert r>=memory+2
 lookup=np.full(K*(1<<n)*8,-1,dtype=np.int32);lookup[s.keys]=np.arange(len(s.keys),dtype=np.int32)
 rr=array('i');cc=array('i');vv=[] if exact else array('d');descriptors=[];cache={}
 def add(i,j,v):assert i>=0;rr.append(i);cc.append(j);vv.append(v if exact else float(v))
 if include_heads:
  cutmap=np.full(K*(1<<(n-1)),-1,dtype=np.int32)
  for i,pk in enumerate(s.keys):
   pk=int(pk);label=pk&7
   if label>1:continue
   key=pk>>3;lev=key>>n;w=key&((1<<n)-1);first=1 if w>>(n-1)&1 else -1
   cut=(lev<<(n-1))|(w>>1) if label==0 else ((lev+first)<<(n-1))|(w&((1<<(n-1))-1))
   col=int(cutmap[cut])
   if col<0:col=len(descriptors);cutmap[cut]=col;descriptors.append(('H',cut))
   add(i,col,Q(1) if exact else 1.)
  del cutmap
 m=2*r-3
 for anchor,w in windows(K,m):
  lev=anchor>>m;h=m//2
  for side in ('L','R'):
   split='S1' if w[h]==1 else ('S2' if side=='L' else 'S3');labels=(split,side,'M','R' if side=='L' else 'L');states=[];word=w;head=h;valid=True
   for label in labels:
    states.append((word,head,label))
    if label=='L':head-=1
    elif label=='R':head+=1
    else:word,head=literal.literal_move(word,head,label)
    if not literal.admissible(lev,word,K):valid=False;break
   if not valid:continue
   assert (word,head)==(w,h)
   begin=h-memory;local=w[begin:begin+2*memory+1];start=lev+sum(w[:begin]);cachekey=(start,code(local),split)
   if cachekey not in cache:
    trip={'S1':(1,-1,1),'S2':(1,-1,-1),'S3':(-1,-1,1)}[split];target=local[:memory]+trip+local[memory+1:];cache[cachekey]=getattr(s.law,'exact',s.law).cylinder((start,local))/getattr(s.law,'exact',s.law).cylinder((start,target))
   ratio=cache[cachekey];column=len(descriptors);descriptors.append((side,anchor))
   for k,(word,head,label) in enumerate(states):
    left=r-head;right=r-(len(word)-1-head);assert left>=0 and right>=0
    value=ratio if k in (1,2) else Q(1)
    for exterior in itertools.product((1,-1),repeat=left+right):
     prefix=exterior[:left];suffix=exterior[left:];full=prefix+word+suffix;start=lev-sum(prefix)
     if literal.admissible(start,full,K):add(int(lookup[(((start<<n)|code(full))<<3)|LABELS.index(label)]),column,value)
 if exact:return list(zip(rr,cc,vv)),descriptors
 B=csr_matrix((np.asarray(vv),(np.asarray(rr),np.asarray(cc))),shape=(len(s.keys),len(descriptors)));B.sum_duplicates();return B,descriptors
