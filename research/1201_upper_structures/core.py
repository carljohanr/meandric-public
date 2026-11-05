"""Generalized sparse row builder (copy/extension of 1007 gate_target.py).
Imports the old evaluator (head.Evaluator/lift) and literal moves read-only; all outputs stay in this folder."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ.setdefault(name,'2')
import importlib.util,json,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
HERE=Path(__file__).resolve().parent;RES=HERE.parent
def _load(name,rel):
    spec=importlib.util.spec_from_file_location(name,RES/rel);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
head=_load('head','1005_head_edge_target/search.py')
literal=_load('literal','1005_head_edge_independent_audit/literal.py')

def edit_fields(d):
    width=2*d+1;mask=(1<<width)-1;mapping={};count=0
    for w in range(mask+1):
        for tri in ((1,2,4) if not ((w>>d)&1) else (3,5,6)):
            key=8*w+tri
            if key in mapping:continue
            wr=int(f'{w:0{width}b}'[::-1],2);tr=int(f'{tri:03b}'[::-1],2)
            for a,b in ((w,tri),(mask-w,7-tri),(wr,tr),(mask-wr,7-tr)):mapping[8*a+b]=count
            count+=1
    return mapping,count

KIND_L,KIND_R,KIND_SPLIT,KIND_MERGE=0,1,2,3

class Model:
    """Variables: [U orbits | F (2^k) | E head cuts | G edit | extra...].  extra: object with n (count), terms(w,h0,t,h,kind)->dict, init(), scale"""
    def __init__(self,r,k,q,d,extra=None,meta=False):
        self.r,self.k,self.q,self.d=r,k,q,d
        t0=time.monotonic()
        self.ev=ev=head.Evaluator(r,min(k,r),q);self.emap,self.nd=edit_fields(d) if d>0 else ({},0)
        self.nbase=ev.nu+(1<<k)+ev.ne;self.extra=extra;self.nx=extra.n if extra else 0
        self.nv=self.nbase+self.nd+self.nx
        fullgroups=1<<(2*r+2);h0=r+1
        words=np.arange(fullgroups,dtype=np.uint64);partners=head.reverse(words,2*r+3);partners=np.minimum(partners,2*fullgroups-1-partners)
        reps=words[partners>=words];self.multiplicity=np.where(partners[reps]==reps,1.,2.);self.reps=reps;ng=self.ng=len(reps)
        cap=(5+(extra.maxnnz if extra else 0))
        vals=np.empty(cap*8*ng,dtype=np.float64);cols=np.empty(cap*8*ng,dtype=np.int32)
        ptr=np.empty(8*ng+1,dtype=np.int64);own=np.empty(8*ng,dtype=np.int32);kinds=np.empty(8*ng,dtype=np.int8) if meta else None
        ptr[0]=0;term=nnz=0;nu=ev.nu;kk=k
        def feature(w,h):
            a=Counter({int(ev.lookup[int(w[h-r:h+r+1],2)]):1})
            for side in (w[:h],w[h+1:][::-1]):
                for j in range(len(side)-kk+1):a[nu+int(side[j:j+kk],2)]+=1
            return a
        for group,g in enumerate(reps):
            physical=f'{g:0{2*r+3}b}'
            for w in (physical,physical.translate(str.maketrans('01','10'))):
                src=feature(w,h0)
                for t,h in literal.outgoing((w,h0)):
                    a=feature(t,h)
                    for j,v in src.items():a[j]-=v
                    if t==w:
                        right=h>h0;kind=KIND_R if right else KIND_L
                        if q:
                            begin=h0-q+right;cut=int(physical[begin:begin+2*q],2)
                            j,sign=ev.field[str(cut)]
                            if j>=0:a[nu+(1<<kk)+j]+=sign if right else -sign
                    else:
                        split=len(t)>len(w);kind=KIND_SPLIT if split else KIND_MERGE
                        if d:
                            short,sh=(w,h0) if split else (t,h);long,lh=(t,h) if split else (w,h0)
                            window=int(short[sh-d:sh+d+1],2);tri=int(long[lh-1:lh+2],2)
                            a[self.nbase+self.emap[8*window+tri]]+=1 if split else -1
                    if extra:
                        for j,v in extra.terms(w,h0,t,h,kind,physical).items():a[self.nbase+self.nd+j]+=v
                    for j,v in a.items():
                        if v:
                            vals[nnz]=v;cols[nnz]=j;nnz+=1
                    own[term]=group;
                    if meta:kinds[term]=kind
                    term+=1;ptr[term]=nnz
        self.X=csr_matrix((vals[:nnz],cols[:nnz],ptr[:term+1]),shape=(term,self.nv),copy=False);self.owner=own[:term].copy()
        self.kinds=kinds[:term].copy() if meta else None
        self.setup_seconds=time.monotonic()-t0
        nu_,ne=ev.nu,ev.ne
        self.scale=np.r_[np.full(nu_,np.sqrt(nu_/(1<<k))),np.ones(1<<k),np.full(ne,np.sqrt(max(1,ne)/(1<<k))),np.full(self.nd,np.sqrt(max(1,self.nd)/(1<<k))),(extra.scale if extra else np.zeros(0))]
        lo=np.r_[[-24]*nu_,[-3]*(1<<k),[-3]*(ne+self.nd),(extra.lo if extra else [])];hi=np.r_[[12]*nu_,[2]*(1<<k),[3]*(ne+self.nd),(extra.hi if extra else [])]
        self.lo,self.hi=np.asarray(lo,float),np.asarray(hi,float)
    def values(self,x):
        ex=np.exp(self.X@x);sums=np.bincount(self.owner,weights=ex,minlength=self.ng)/2;return ex,sums
    def x_from_witness(self,old):
        """old: json of a head/edit witness with r<=self.r etc. Edit/extra variables from G if same d else zero."""
        base=head.lift(old,self.ev);x=np.r_[base,np.zeros(self.nd+self.nx)]
        if self.d and old.get('d')==self.d:
            for key,j in self.emap.items():x[self.nbase+j]=np.log(float(Q(old['G'][old['edit_lookup'][str(key)]])))
        return x
