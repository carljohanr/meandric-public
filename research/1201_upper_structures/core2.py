"""Model2: like core.Model but row half-length `a` is independent of the U radius rU.
Row = word of length 2a+1 around the head (head at index a).  Old model: a=r+1, rU=r.
Families: U (window radius rU, reflection orbits), F (k-grams, k<=a), E (head cut half-width q<=a), G (edit, d<=a-1),
optional refinement V (window radius rU+1 only for windows whose inner radius-rU window is in `refine_inner`),
optional generic `extra` hook."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ.setdefault(name,'1')
import time
from collections import Counter
from core import head,literal,edit_fields,KIND_L,KIND_R,KIND_SPLIT,KIND_MERGE,np,csr_matrix,Q,json

def rev_int(v,width):return int(f'{v:0{width}b}'[::-1],2)

class Model2:
    def __init__(self,a,rU,k,q,d,refine=None,meta=False,extra=None):
        assert rU<=a-1 and k<=a and q<=a and d<=a-1
        self.a,self.rU,self.k,self.q,self.d=a,rU,k,q,d
        t0=time.monotonic()
        reps,lookup=head.structure(rU)         # window radius rU orbit structure
        self.nu=len(reps);self.lookup=lookup
        self.field,_,self.ne=head.fields(q)
        self.emap,self.nd=edit_fields(d) if d>0 else ({},0)
        self.refine=refine or {}                 # window_{rU+1} int -> var idx (relative to refine block)
        self.nr=(max(self.refine.values())+1) if self.refine else 0
        self.inner_set=set(w>>1&((1<<(2*rU+1))-1) for w in self.refine) if self.refine else set()
        self.oU=0;self.oF=self.nu;self.oE=self.oF+(1<<k);self.oG=self.oE+self.ne;self.oR=self.oG+self.nd;self.nv=self.oR+self.nr
        self.nx=extra.n if extra else 0;self.oX=self.nv;self.nv+=self.nx
        L=2*a+1
        fullgroups=1<<(L-1)
        words=np.arange(fullgroups,dtype=np.uint64);partners=head.reverse(words,L)
        partners=np.minimum(partners,2*fullgroups-1-partners)
        sel=partners>=words;reps_rows=words[sel];self.multiplicity=np.where(partners[reps_rows]==reps_rows,1.,2.);self.rows=reps_rows;ng=self.ng=len(reps_rows)
        cap=8
        vals=np.empty(cap*8*ng,dtype=np.float64);cols=np.empty(cap*8*ng,dtype=np.int32)
        ptr=np.empty(8*ng+1,dtype=np.int64);own=np.empty(8*ng,dtype=np.int32);kinds=np.empty(8*ng,dtype=np.int8)
        ptr[0]=0;term=nnz=0
        rU_=rU;kk=k;oF=self.oF;oE=self.oE;oG=self.oG;oR=self.oR;refine=self.refine;lookup_=lookup;rV=rU+1
        def feature(w,h):
            idx=int(lookup_[int(w[h-rU_:h+rU_+1],2)])
            c=Counter({idx:1})
            if refine:
                key=int(w[h-rV:h+rV+1],2) if h-rV>=0 and h+rV+1<=len(w) else -1
                if key in refine:c[oR+refine[key]]+=1
            for side in (w[:h],w[h+1:][::-1]):
                for j in range(len(side)-kk+1):c[oF+int(side[j:j+kk],2)]+=1
            return c
        for group,g in enumerate(reps_rows):
            physical=f'{g:0{L}b}'
            for w in (physical,physical.translate(str.maketrans('01','10'))):
                src=feature(w,a)
                for t,h in literal.outgoing((w,a)):
                    ft=feature(t,h)
                    for j,v in src.items():ft[j]-=v
                    if t==w:
                        right=h>a;kind=KIND_R if right else KIND_L
                        if q:
                            begin=a-q+right;cut=int(physical[begin:begin+2*q],2)
                            j,sign=self.field[str(cut)]
                            if j>=0:ft[oE+j]+=sign if right else -sign
                    else:
                        split=len(t)>len(w);kind=KIND_SPLIT if split else KIND_MERGE
                        if d:
                            short,sh=(w,a) if split else (t,h);long,lh=(t,h) if split else (w,a)
                            window=int(short[sh-d:sh+d+1],2);tri=int(long[lh-1:lh+2],2)
                            ft[oG+self.emap[8*window+tri]]+=1 if split else -1
                    if extra:
                        for j,v in extra.terms(w,a,t,h,kind,physical).items():ft[self.oX+j]+=v
                    for j,v in ft.items():
                        if v:
                            assert nnz<len(vals);vals[nnz]=v;cols[nnz]=j;nnz+=1
                    own[term]=group;kinds[term]=kind;term+=1;ptr[term]=nnz
        self.X=csr_matrix((vals[:nnz],cols[:nnz],ptr[:term+1]),shape=(term,self.nv),copy=False);self.owner=own[:term].copy();self.kinds=kinds[:term].copy()
        self.setup_seconds=time.monotonic()-t0
    def values(self,x):
        ex=np.exp(self.X@x);sums=np.bincount(self.owner,weights=ex,minlength=self.ng)/2;return ex,sums

def refine_map(rU,inner_windows):
    """variables for window_{rU+1} = (a,inner,b), inner in inner_windows (ints of width 2rU+1); reflection ties (a,inner,b)~(b,rev(inner),a)."""
    wI=2*rU+1;mp={};n=0
    for inner in sorted(inner_windows):
        for x in (0,1):
            for y in (0,1):
                key=(x<<(wI+1))|(inner<<1)|y
                if key in mp:continue
                rk=(y<<(wI+1))|(rev_int(inner,wI)<<1)|x
                mp[key]=n;mp[rk]=n;n+=1
    return mp
