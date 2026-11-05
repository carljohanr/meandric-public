"""Explicit local stationary directions from actual two- and four-move cycles.

Head cycles use a shared cut intensity. Edit cycles are S,L,M,R or S,R,M,L
on a core of length 2r-3. Every state of such a cycle fits in the radius-r
window. Give each edge core density beta(anchor)/beta(current core); summing
over exterior completions makes its physical mass exactly beta(anchor).
Thus incoming/outgoing flows cancel and the two count drifts sum to zero.
The implementation must pass an exact original-row control before use.
"""
import sys,itertools
from pathlib import Path
from array import array
from types import SimpleNamespace
import numpy as np
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parent.parent
sys.path[:0]=[str(ROOT/'971_k14_streaming_preflight'),str(ROOT/'945_lower_scaling')]
from stream_builder import windows,decode,code,LABELS
from float_system import FloatLaw
import importlib.util
spec=importlib.util.spec_from_file_location('literal_cycle',ROOT/'940_integer_literal_verifier/verify_sorted.py');literal=importlib.util.module_from_spec(spec);spec.loader.exec_module(literal)

def geometry(c,r):
    K=c['K'];oldn=2*c['r']+1;n=2*r+1;extra=r-c['r'];assert extra>=0
    law=FloatLaw(K,c['up'],c['kappa']);tab=dict(zip(c['flow_keys'],c['j']))
    keys=array('I');vw=array('i');bw=array('d');base=array('d');da=array('b');db=array('b')
    for wi,(wk,w) in enumerate(windows(K,n)):
        lev=wk>>n;bw.append(law.cylinder((lev,w)));oldword=w[extra:n-extra] if extra else w;oldlev=lev+sum(w[:extra]);oldcode=code(oldword)
        for label in literal.legal_labels(lev,w,r,K):
            keys.append((wk<<3)|LABELS.index(label));vw.append(wi);base.append(tab[oldlev,oldcode,label]);d=literal.DRIFT[label];da.append(d[0]);db.append(d[1])
    return SimpleNamespace(K=K,r=r,keys=np.asarray(keys),vw=np.asarray(vw),bw=np.asarray(bw),base=np.asarray(base),da=np.asarray(da),db=np.asarray(db),law=law)

def kernel(s,exact=False,include_heads=True):
    K,r=s.K,s.r;n=2*r+1;law=s.law.exact if exact else s.law
    lookup={int(k):i for i,k in enumerate(s.keys)};rr=array('i');cc=array('i');vv=[] if exact else array('d');descriptors=[]
    def add(i,j,v):rr.append(i);cc.append(j);vv.append(v if exact else float(v))
    if include_heads:
        cutmap={}
        for i,pk in enumerate(s.keys):
            pk=int(pk);label=pk&7
            if label not in (0,1):continue
            key=pk>>3;lev=key>>n;w=key&((1<<n)-1);first=1 if w>>(n-1)&1 else -1
            cut=(lev<<(n-1))|(w>>1) if label==0 else ((lev+first)<<(n-1))|(w&((1<<(n-1))-1))
            if cut not in cutmap:cutmap[cut]=len(descriptors);descriptors.append(('H',cut))
            add(i,cutmap[cut],1)
    m=2*r-3
    for anchor,w in windows(K,m):
        lev=anchor>>m;h=m//2;mass=law.cylinder((lev,w))
        for side in ('L','R'):
            split='S1' if w[h]==1 else ('S2' if side=='L' else 'S3')
            labels=(split,side,'M','R' if side=='L' else 'L');states=[];word=w;head=h;valid=True
            for label in labels:
                states.append((word,head,label))
                if label=='L':word,head=word,head-1
                elif label=='R':word,head=word,head+1
                else:word,head=literal.literal_move(word,head,label)
                if not literal.admissible(lev,word,K):valid=False;break
            if not valid:continue
            assert (word,head)==(w,h)
            column=len(descriptors);descriptors.append((side,anchor))
            for word,head,label in states:
                left=r-head;right=r-(len(word)-1-head);assert left>=0 and right>=0
                ratio=mass/law.cylinder((lev,word))
                for exterior in itertools.product((1,-1),repeat=left+right):
                    prefix=exterior[:left];suffix=exterior[left:];full=prefix+word+suffix;start=lev-sum(prefix)
                    if literal.admissible(start,full,K):
                        pk=(((start<<n)|code(full))<<3)|LABELS.index(label)
                        add(lookup[pk],column,ratio)
    if exact:return list(zip(rr,cc,vv)),descriptors
    B=csr_matrix((np.asarray(vv),(np.asarray(rr),np.asarray(cc))),shape=(len(s.keys),len(descriptors)));B.sum_duplicates()
    return B,descriptors
