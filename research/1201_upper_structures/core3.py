"""Model3: base families (rU,k,q,d) at row half-length a, plus selected refinement of each family to the next level
(rU+1,k+1,q+1,d+1) as multiplicative correction variables, only for selected parents.  a >= base+1 + 1 when refinement used."""
import os
import time
from collections import Counter
from core import head,literal,edit_fields,KIND_L,KIND_R,KIND_SPLIT,KIND_MERGE,np,csr_matrix,Q,json
from core2 import rev_int

class Model3:
    def __init__(self,a,rU,k,q,d,selU=None,selE=None,selG=None,selF=None,meta=False,cap=None):
        """selU: set of inner window ints (width 2rU+1) ; selE: set of parent cut ints (width 2q) ; selG: set of (window_d,tri) keys 8*w+tri ; selF: set of parent gram ints (width k)."""
        refined=any(s for s in (selU,selE,selG,selF))
        self.a,self.rU,self.k,self.q,self.d=a,rU,k,q,d
        assert a>=max(rU+1+(selU is not None and len(selU)>0),k+(1 if selF else 0),q+(1 if selE else 0),d+1+(1 if selG else 0))
        t0=time.monotonic()
        reps,lookup=head.structure(rU);self.nu=len(reps);self.lookup=lookup
        field,_,ne=head.fields(q) if q else ({},None,0)
        self.field=field;self.ne=ne;emap,nd=edit_fields(d) if d else ({},0)
        self.emap=emap;self.nd=nd
        # fine tables
        self.selU=selU or set();self.selE=selE or set();self.selG=selG or set();self.selF=selF or set()
        # U fine
        from core2 import refine_map
        self.rmU=refine_map(rU,self.selU) if self.selU else {};self.nrU=(max(self.rmU.values())+1) if self.rmU else 0
        # E fine: orbits of fine cuts whose parent in selE
        self.fE={};self.nrE=0
        if self.selE:
            fieldF,_,_=head.fields(q+1);w=2*(q+1);m=(1<<w)-1
            for cut in range(1<<w):
                par=(cut>>1)&((1<<(2*q))-1)
                if par not in self.selE:continue
                j,s=fieldF[str(cut)]
                if j<0:continue
                self.fE[cut]=(j,s)
            # remap j to compact
            ids=sorted(set(j for j,s in self.fE.values()));mp={j:i for i,j in enumerate(ids)}
            self.fE={c:(mp[j],s) for c,(j,s) in self.fE.items()};self.nrE=len(ids)
        # G fine
        self.fG={};self.nrG=0
        if self.selG:
            emapF,ndF=edit_fields(d+1);ids={}
            for key,j in emapF.items():
                w,tri=divmod(key,8);par=(w>>1)&((1<<(2*d+1))-1)
                if 8*par+tri not in self.selG:continue
                if j not in ids:ids[j]=len(ids)
                self.fG[key]=ids[j]
            self.nrG=len(ids)
        # F fine: grams of length k+1 whose parent (last k letters, nearest head) in selF; var per gram
        self.fF={};self.nrF=0
        if self.selF:
            for g in range(1<<(k+1)):
                par=g&((1<<k)-1)      # last k letters (nearest head) = low bits
                if par in self.selF:self.fF[g]=len(self.fF)
            self.nrF=len(self.fF)
        o=self.nu;self.oF=o;o+=1<<k;self.oE=o;o+=ne;self.oG=o;o+=nd
        self.oRU=o;o+=self.nrU;self.oRE=o;o+=self.nrE;self.oRG=o;o+=self.nrG;self.oRF=o;o+=self.nrF;self.nv=o
        L=2*a+1;fullgroups=1<<(L-1)
        words=np.arange(fullgroups,dtype=np.uint64);partners=head.reverse(words,L);partners=np.minimum(partners,2*fullgroups-1-partners)
        sel=partners>=words;rows=words[sel];self.multiplicity=np.where(partners[rows]==rows,1.,2.);self.rows=rows;ng=self.ng=len(rows)
        cap=cap or (12 if refined else 5);vals=np.empty(8*cap*ng,dtype=np.float64);cols=np.empty(8*cap*ng,dtype=np.int32)
        ptr=np.empty(8*ng+1,dtype=np.int64);own=np.empty(8*ng,dtype=np.int32);kinds=np.empty(8*ng,dtype=np.int8)
        ptr[0]=0;term=nnz=0;rV=rU+1
        lookup_=lookup;kk=k
        def feature(w,h):
            idx=int(lookup_[int(w[h-rU:h+rU+1],2)]);c=Counter({idx:1})
            if self.rmU:
                key=int(w[h-rV:h+rV+1],2) if h-rV>=0 and h+rV+1<=len(w) else -1
                if key in self.rmU:c[self.oRU+self.rmU[key]]+=1
            for side in (w[:h],w[h+1:][::-1]):
                for j in range(len(side)-kk+1):c[self.oF+int(side[j:j+kk],2)]+=1
                if self.fF:
                    for j in range(len(side)-kk):
                        g=int(side[j:j+kk+1],2)
                        # side string runs far->near ; last letters are nearest head ; parent = last k letters of the (k+1)-gram
                        if g in self.fF:c[self.oRF+self.fF[g]]+=1
            return c
        for group,g in enumerate(rows):
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
                            j,sign=field[str(cut)]
                            if j>=0:ft[self.oE+j]+=sign if right else -sign
                        if self.fE:
                            begin=a-(q+1)+right;cut=int(physical[begin:begin+2*q+2],2)
                            if cut in self.fE:
                                j,sign=self.fE[cut];ft[self.oRE+j]+=sign if right else -sign
                    else:
                        split=len(t)>len(w);kind=KIND_SPLIT if split else KIND_MERGE
                        short,sh=(w,a) if split else (t,h);long,lh=(t,h) if split else (w,a)
                        tri=int(long[lh-1:lh+2],2)
                        if d:
                            window=int(short[sh-d:sh+d+1],2)
                            ft[self.oG+emap[8*window+tri]]+=1 if split else -1
                        if self.fG:
                            window=int(short[sh-d-1:sh+d+2],2)
                            key=8*window+tri
                            if key in self.fG:ft[self.oRG+self.fG[key]]+=1 if split else -1
                    for j,v in ft.items():
                        if v:
                            vals[nnz]=v;cols[nnz]=j;nnz+=1
                    own[term]=group;kinds[term]=kind;term+=1;ptr[term]=nnz
        self.X=csr_matrix((vals[:nnz],cols[:nnz],ptr[:term+1]),shape=(term,self.nv),copy=False);self.owner=own[:term].copy();self.kinds=kinds[:term].copy()
        self.setup_seconds=time.monotonic()-t0
    def values(self,x):
        ex=np.exp(self.X@x);sums=np.bincount(self.owner,weights=ex,minlength=self.ng)/2;return ex,sums
