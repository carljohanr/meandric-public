"""Original density rows from integer move topology and exact local ratios.

All exterior extension probabilities cancel against the target cylinder.
Incoming coefficient is1 for heads, beta(before)/beta(after) for a local
split/merge, independent of the exterior completion. Ratios are computed
exactly on the central2m+1 or2m+3 edges and converted once.
"""
import sys,time
from pathlib import Path
from array import array
from types import SimpleNamespace
from fractions import Fraction as Q
import numpy as np
from scipy.sparse import csr_matrix
from memory_two import FloatMemoryTwo,ROOT
from memory_order import FloatOrderLaw
sys.path.insert(0,str(ROOT/'971_k14_streaming_preflight'))
from stream_builder import windows,decode,code,LABELS,MIRROR,reflection_key
from float_system import FloatLaw
from band_model import strip_ok

def build_local(K,r,up,kappa,parameters=(),exact_control=False,law_override=None,memory_override=None,retain_topology=False):
 start=time.monotonic();parameters=tuple(parameters);m=memory_override if memory_override is not None else len(parameters)+1;assert r>=m+1
 law=FloatLaw(K,up,kappa) if not parameters else (FloatMemoryTwo(K,up,kappa,parameters[0]) if len(parameters)==1 else FloatOrderLaw(K,up,kappa,parameters));law=law if law_override is None else law_override;E=law.exact;n=2*r+1
 win_keys=array('I');bw=array('d')
 for key,w in windows(K,n):win_keys.append(key);bw.append(float(law.cylinder((key>>n,w))))
 wi={int(key):i for i,key in enumerate(win_keys)};row_keys=array('I',[key for key,_ in windows(K,n+2)]);ri={int(key):i for i,key in enumerate(row_keys)};inner=np.empty(len(row_keys),dtype=np.int32)
 for row,key in enumerate(row_keys):
  lev=int(key)>>(n+2);word=int(key)&((1<<(n+2))-1);first=1 if word>>(n+1)&1 else -1;inner[row]=wi[((lev+first)<<n)|((word>>1)&((1<<n)-1))]
 outgoing_order=np.argsort(inner,kind='stable');out_offsets=np.r_[0,np.cumsum(np.bincount(inner,minlength=len(win_keys)))];cache={};ratio_ids={};ratio_keys=[None];ratio_type=array('i');incoming_flags=array('b');outgoing_flags=array('b');rr=array('i');cc=array('i');vv=array('d');keys=array('I');vw=array('i');da=array('b');db=array('b');exactrows=[{} for _ in row_keys] if exact_control else None
 def add(row,col,value,pflag=0,nflag=0):
  if not value and not retain_topology:return
  rr.append(row);cc.append(col);vv.append(float(value))
  if retain_topology:incoming_flags.append(pflag);outgoing_flags.append(nflag)
  if exactrows is not None:exactrows[row][col]=value
 def coefficient(lev,w,label):
  if label<2:return Q(1),0
  begin=r-m if label<5 else r-m-1;length=2*m+1 if label<5 else 2*m+3;local=w[begin:begin+length];a=lev+sum(w[:begin]);key=(a,code(local),label)
  if key not in cache:
   if label<5:
    trip=((1,-1,1),(1,-1,-1),(-1,-1,1))[label-2];target=local[:m]+trip+local[m+1:]
   else:target=local[:m]+(local[m]+local[m+2]+1,)+local[m+3:]
   cache[key]=E.cylinder((a,local))/E.cylinder((a,target));ratio_ids[key]=len(ratio_keys);ratio_keys.append(key)
  return cache[key],ratio_ids[key]
 for i,key in enumerate(win_keys):
  key=int(key);lev=key>>n;w=decode(key&((1<<n)-1),n);end=lev+sum(w);left=[(f2,f1) for f1 in (1,-1) if 0<=lev-f1<K for f2 in (1,-1) if 0<=lev-f1-f2<K];right=[(g1,g2) for g1 in (1,-1) if 0<=end+g1<K for g2 in (1,-1) if 0<=end+g1+g2<K]
  moves=[(0,(-1,1),[(lev-f1-f2,(f2,f1)+w) for f2,f1 in left]),(1,(1,-1),[(lev,w+(g1,g2)) for g1,g2 in right])]
  for label,trip in ([(2,(1,-1,1))] if w[r]==1 else [(3,(1,-1,-1)),(4,(-1,-1,1))]):
   target=w[:r]+trip+w[r+1:]
   if strip_ok(lev,target,K):moves.append((label,(1,1),[(lev,target)]))
  if w[r]==1 and (w[r-1],w[r+1])!=(1,1):
   core=w[:r-1]+(w[r-1]+w[r+1]+1,)+w[r+2:];assert core[:m]==w[:m] and core[-m:]==w[-m:];moves.append((5,(-1,-1),[(lev-f1-f2,(f2,f1)+core+(g1,g2)) for f2,f1 in left for g1,g2 in right]))
  for label,d,targets in moves:
   col=len(keys);keys.append((key<<3)|label);vw.append(i);da.append(d[0]);db.append(d[1]);c,rid=coefficient(lev,w,label);ratio_type.append(rid);incoming={ri[(a<<(n+2))|code(word)] for a,word in targets}
   for row in incoming:add(row,col,c-1 if inner[row]==i else c,1,int(inner[row]==i))
   for row in outgoing_order[out_offsets[i]:out_offsets[i+1]]:
    if row not in incoming:add(int(row),col,Q(-1),0,1)
 M=csr_matrix((np.frombuffer(vv,dtype=np.float64),(np.frombuffer(rr,dtype=np.int32),np.frombuffer(cc,dtype=np.int32))),shape=(len(row_keys),len(keys)));topology=None
 if retain_topology:
  P=csr_matrix((np.asarray(incoming_flags,dtype=np.int8),(np.asarray(rr),np.asarray(cc))),shape=M.shape);N=csr_matrix((np.asarray(outgoing_flags,dtype=np.int8),(np.asarray(rr),np.asarray(cc))),shape=M.shape);assert np.array_equal(P.indices,M.indices) and np.array_equal(N.indices,M.indices)
  topology=(P.data.copy(),N.data.copy(),np.asarray(ratio_type),ratio_keys)
 flow_index={int(key):i for i,key in enumerate(keys)};flowperm=np.array([flow_index[(reflection_key(int(key)>>3,n,K)<<3)|MIRROR[int(key)&7]] for key in keys],dtype=np.int32);rowperm=np.array([ri[reflection_key(int(key),n+2,K)] for key in row_keys],dtype=np.int32)
 return SimpleNamespace(K=K,r=r,M=M,keys=np.asarray(keys),vw=np.asarray(vw),bw=np.asarray(bw),da=np.asarray(da),db=np.asarray(db),flowperm=flowperm,rowperm=rowperm,row_keys=np.asarray(row_keys),exactrows=exactrows,law=law,local_ratio_cache=cache,topology=topology,memory=m,seconds=time.monotonic()-start)
