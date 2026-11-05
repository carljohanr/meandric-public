"""Literal deep cycles and exact sparse entropy increments without a full fine model."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
from fractions import Fraction as Q
from array import array
import json,itertools,heapq,time
import numpy as np
from local_cycle_family import literal,windows,code,LABELS,FloatOrderLaw
from scipy.sparse import coo_matrix

def phi(x):
 out=(1+x)*np.log1p(x)-x;small=abs(x)<1e-4;t=x[small];out[small]=t*t*(.5+t*(-1/6+t*(1/12+t*(-1/20+t/30))));return out

class SparseFineBaseline:
 def __init__(self,source,r):
  self.source=Path(source);self.r=r;self.n=2*r+1
  import hashlib
  cursor=Path(source);overlays=[];top=None
  while True:
   with np.load(cursor) as z:
    meta=json.loads(str(z['metadata']))
    if 'physical' in z:self.metadata=meta;keys=z['keys'].astype(np.int64);F=z['physical'].copy();break
    assert meta['status']=='NUMERICAL_ADAPTIVE_SPARSE_DEEP_CYCLE_NOT_CERTIFIED' and r>=int(meta['r']) and meta['minimum_changed_density']>0
    kk=z['changed_keys'].astype(np.int64);dd=z['density_increment'].copy();assert np.all(kk[1:]>kk[:-1]) and np.isfinite(dd).all();ww,ii=np.unique(kk>>3,return_inverse=True);overlays.append(dict(r=int(meta['r']),keys=kk,delta=dd,windows=ww,total=np.bincount(ii,weights=dd)))
    if top is None:top=meta
   parent=Path(meta['source'])
   if not parent.exists():parent=Path(__file__).resolve().parent.parent.parent/parent
   assert hashlib.file_digest(parent.open('rb'),'sha256').hexdigest()==meta['source_sha256'];cursor=parent
  self.overlays=list(reversed(overlays))
  meta=self.metadata;self.K=int(meta['K']);self.r0=int(meta['r']);self.n0=2*self.r0+1;assert r>=self.r0
  up={int(i):Q(v) for i,v in meta['up'].items()};params=tuple(map(Q,meta['memory_parameters']));self.parameters=params;self.law=FloatOrderLaw(self.K,up,Q(meta['kappa']),params)
  if 'memory_profile' in meta:
   from mixed_profile_law import make_law
   self.law=make_law(self.K,up,Q(meta['kappa']),params,{int(g):Q(v) for g,v in meta['memory_profile'].items()})
  self.m=self.law.memory
  self.state=np.zeros(self.K*(1<<self.m));self.forward=np.zeros(self.K*(1<<(self.m+1)))
  def encode(v):return (v[0]<<(len(v)-1))+sum(int(b>a)<<(len(v)-2-i) for i,(a,b) in enumerate(zip(v,v[1:])))
  for v,value in self.law.state.items():self.state[encode(v)]=value
  for v,value in self.law.forward.items():self.forward[encode(v)]=value
  wk=keys>>3;unique,inv=np.unique(wk,return_inverse=True);weights=self.weights(unique,self.n0);density=F/weights[inv];assert density.min()>0
  self.lookup=np.full(self.K*(1<<self.n0)*8,-1,dtype=np.int32);self.lookup[keys]=np.arange(len(keys),dtype=np.int32);self.density=density
  self.total=np.zeros(self.K*(1<<self.n0));self.total[unique]=np.bincount(inv,weights=density)
  self.mass=float(F.sum());self.h=float(-F@np.log(density/self.total[wk])/self.mass)
  if 'mu' in meta:assert abs(np.exp(2*self.h)/float(meta['mu'])-1)<1e-12
  if top is not None:self.mass=float(top['mass']);self.h=float(top.get('entropy',np.log(top['mu'])/2))
  self.initial_mu=float(np.exp(2*self.h))
 def weights(self,wk,n=None):
  n=self.n if n is None else n;wk=np.asarray(wk,dtype=np.int64);lev=wk>>n;word=wk&((1<<n)-1);value=self.state[(lev<<self.m)|(word>>(n-self.m))].copy();height=lev.copy()
  for i in range(self.m,n):
   value*=self.forward[(height<<(self.m+1))|((word>>(n-1-i))&((1<<(self.m+1))-1))];height+=2*((word>>(n-1-(i-self.m)))&1)-1
  return value
 def restrict(self,wk,radius=None):
  wk=np.asarray(wk,dtype=np.int64);word=wk&((1<<self.n)-1);lev=wk>>self.n;radius=self.r0 if radius is None else radius;n0=2*radius+1;extra=self.r-radius;assert extra>=0
  for i in range(extra):lev+=2*((word>>(self.n-1-i))&1)-1
  return (lev<<n0)|((word>>extra)&((1<<n0)-1))
 def add_overlay(self,keys,authority,values,out):
  indices=np.searchsorted(authority,keys);valid=indices<len(authority);valid[valid]=authority[indices[valid]]==keys[valid];out[valid]+=values[indices[valid]]
 def window_total(self,wk):
  total=self.total[self.restrict(wk)].copy()
  for overlay in self.overlays:self.add_overlay(self.restrict(wk,overlay['r']),overlay['windows'],overlay['total'],total)
  return total
 def values(self,keys):
  keys=np.asarray(keys,dtype=np.int64);wk=keys>>3;old=self.restrict(wk);indices=self.lookup[(old<<3)|(keys&7)];assert indices.min()>=0
  density=self.density[indices].copy()
  for overlay in self.overlays:self.add_overlay((self.restrict(wk,overlay['r'])<<3)|(keys&7),overlay['keys'],overlay['delta'],density)
  assert density.min()>0
  return self.weights(wk),density,self.window_total(wk)

class LiteralCycleGenerator:
 def __init__(self,K,r,law,parameters):self.K=K;self.r=r;self.n=2*r+1;self.m=2*r-3;self.law=law;self.memory=len(parameters)+1;self.ratios={};assert r>=self.memory+2
 def descriptors(self,heads=True,edits=True):
  if heads:
   for cut,_ in windows(self.K,self.n-1):yield (0,cut)
  if edits:
   for anchor,w in windows(self.K,self.m):
    for side in (1,2):
     if self.entries((side,anchor),word=w):yield (side,anchor)
 def entries(self,descriptor,exact=False,word=None):
  kind,key=map(int,descriptor);K,r,n=self.K,self.r,self.n
  if kind==0:
   lev=key>>(n-1);bits=key&((1<<(n-1))-1);end=lev+2*bits.bit_count()-(n-1);out=[]
   for step in (1,-1):
    if 0<=end+step<K:out.append((((key<<1)|int(step==1))<<3, Q(1) if exact else 1.))
    start=lev-step
    if 0<=start<K:out.append(((((start<<n)|(int(step==1)<<(n-1))|bits)<<3)|1,Q(1) if exact else 1.))
   return out
  m=self.m;lev=key>>m;bits=key&((1<<m)-1);w=word if word is not None else tuple(1 if (bits>>j)&1 else -1 for j in range(m-1,-1,-1));h=m//2;side='L' if kind==1 else 'R';split='S1' if w[h]==1 else ('S2' if kind==1 else 'S3');labels=(split,side,'M','R' if kind==1 else 'L');states=[];word=w;head=h
  for label in labels:
   states.append((word,head,label))
   if label=='L':head-=1
   elif label=='R':head+=1
   else:word,head=literal.literal_move(word,head,label)
   if not literal.admissible(lev,word,K):return []
  assert (word,head)==(w,h);begin=h-self.memory;local=w[begin:begin+2*self.memory+1];start=lev+sum(w[:begin]);cachekey=(start,code(local),split)
  if cachekey not in self.ratios:
   trip={'S1':(1,-1,1),'S2':(1,-1,-1),'S3':(-1,-1,1)}[split];target=local[:self.memory]+trip+local[self.memory+1:];authority=getattr(self.law,'exact',self.law);self.ratios[cachekey]=authority.cylinder((start,local))/authority.cylinder((start,target))
  ratio=self.ratios[cachekey];out=[]
  for stage,(word,head,label) in enumerate(states):
   left=r-head;right=r-(len(word)-1-head);value=ratio if stage in (1,2) else Q(1)
   if not exact:value=float(value)
   for exterior in itertools.product((1,-1),repeat=left+right):
    prefix=exterior[:left];full=prefix+word+exterior[left:];start=lev-sum(prefix)
    if literal.admissible(start,full,K):out.append(((((start<<n)|code(full))<<3)|LABELS.index(label),value))
  return out

def build_sparse_problem(base,generator,descriptors):
 rr=array('q');cc=array('i');vv=array('d')
 for col,descriptor in enumerate(descriptors):
  for key,value in generator.entries(descriptor):rr.append(key);cc.append(col);vv.append(value)
 keys,row=np.unique(np.asarray(rr),return_inverse=True);B=coo_matrix((np.asarray(vv),(row,np.asarray(cc))),shape=(len(keys),len(descriptors))).tocsr();B.sum_duplicates();wk,vertex=np.unique(keys>>3,return_inverse=True);beta,j,J=base.values(keys);wbeta=base.weights(wk);wJ=base.window_total(wk);return dict(keys=keys,B=B,vertex=vertex,beta=beta,j=j,nu=wbeta*wJ,logs=np.log(j/J),window_keys=wk)

def sparse_entropy(base,problem,coefficients):
 B=problem['B'];delta=B@coefficients;j=problem['j'];beta=problem['beta'];df=beta*delta;dn=np.bincount(problem['vertex'],weights=df,minlength=len(problem['nu']));mass=base.mass+df.sum();assert mass>0 and np.all(j+delta>0) and np.all(problem['nu']+dn>0)
 gain=(df@(-problem['logs']-base.h)-(beta*j)@phi(delta/j)+problem['nu']@phi(dn/problem['nu']))/mass;h=base.h+gain;logs=problem['logs']+np.log1p(delta/j)-np.log1p(dn/problem['nu'])[problem['vertex']];gradient=np.asarray(B.T@(beta*(-logs-h))).ravel()/mass
 return h,gradient,delta

def sample_window_keys(base,length,count,rng,law_fraction=.7):
 """Mixture of reference-law cylinders and uniform legal height walks."""
 assert length>=base.m and 0<=law_fraction<=1;n_law=int(count*law_fraction);n_uniform=count-n_law;pieces=[]
 if n_law:
  packed=rng.choice(len(base.state),size=n_law,p=base.state/base.state.sum());start=packed>>base.m;word=packed&((1<<base.m)-1);current=start.copy();tail=word.copy()
  for _ in range(base.m,length):
   indices=(current<<(base.m+1))|(tail<<1);up=base.forward[indices|1];down=base.forward[indices];assert np.max(abs(up+down-1))<1e-12
   bit=(rng.random(n_law)<up).astype(np.int64);word=(word<<1)|bit;current+=2*(tail>>(base.m-1))-1;tail=((tail<<1)|bit)&((1<<base.m)-1);assert current.min()>=0 and current.max()<base.K
  pieces.append((start<<length)|word)
 if n_uniform:
  dp=np.ones((length+1,base.K))
  for L in range(1,length+1):dp[L]=np.r_[dp[L-1,1:],0]+np.r_[0,dp[L-1,:-1]]
  start=rng.choice(base.K,size=n_uniform,p=dp[length]/dp[length].sum());height=start.copy();word=np.zeros(n_uniform,dtype=np.int64)
  for i in range(length):
   remaining=length-i-1;up=np.zeros(n_uniform);valid=height+1<base.K;up[valid]=dp[remaining,height[valid]+1]/dp[remaining+1,height[valid]];bit=(rng.random(n_uniform)<up).astype(np.int64);word=(word<<1)|bit;height+=2*bit-1;assert height.min()>=0 and height.max()<base.K
  pieces.append((start<<length)|word)
 out=np.concatenate(pieces);rng.shuffle(out);return out

def child_descriptors(descriptor,r,K):
 """One extra tape letter at each end; literal columns refine the old cylinder."""
 kind,key=map(int,descriptor);length=2*r if kind==0 else 2*r-3;lev=key>>length;bits=key&((1<<length)-1);end=lev+2*bits.bit_count()-length
 for prefix in (1,-1):
  start=lev-prefix
  if not 0<=start<K:continue
  for suffix in (1,-1):
   if 0<=end+suffix<K:yield kind,(start<<(length+2))|(int(prefix==1)<<(length+1))|(bits<<1)|int(suffix==1)
