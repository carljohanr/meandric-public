"""Rational reversible quartet law with unchanged stationary triple marginals."""
import sys
from pathlib import Path
from fractions import Fraction as Q
from types import SimpleNamespace
from array import array
import numpy as np
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parent.parent
sys.path[:0]=[str(ROOT/'971_k14_streaming_preflight'),str(ROOT/'914_level_memory_gate')]
from edge_model import EdgeLaw
from stream_builder import windows,decode,code,LABELS,MIRROR,reflection_key
from band_model import strip_ok

class MemoryTwoLaw:
    def __init__(self,K,up,kappa,eta):
        self.K=K;self.eta=Q(eta);assert abs(self.eta)<1
        base=EdgeLaw(K,up,kappa);self.rho=base.rho;self.edges=base.edges
        self.pi={(a,b,c):base.rho[a,b]*base.P(a,b,c)
            for a,b in base.edges for c in (b-1,b+1) if 0<=c<K}
        quartet={(a,b,c,d):p*base.P(b,c,d) for (a,b,c),p in self.pi.items()
            for d in (c-1,c+1) if 0<=d<K}
        minimum={}
        for (a,b,c,d),p in quartet.items():
            if 0<b<K-1 and 0<c<K-1:
                group=min(b,c,K-1-b,K-1-c)
                minimum[group]=min(minimum.get(group,p),p)
        for key in list(quartet):
            a,b,c,d=key
            if 0<b<K-1 and 0<c<K-1:
                group=min(b,c,K-1-b,K-1-c)
                sign=(1 if a<b else -1)*(1 if d<c else -1)
                quartet[key]+=self.eta*minimum[group]*sign
        for key,p in quartet.items():
            assert p>0 and quartet[key[::-1]]==p
            assert quartet[tuple(K-1-x for x in key)]==p
        for (a,b,c),p in self.pi.items():
            assert sum(quartet[a,b,c,d] for d in (c-1,c+1) if 0<=d<K)==p
            assert sum(quartet[z,a,b,c] for z in (a-1,a+1) if 0<=z<K)==p
        self.quartet=quartet
        self.forward={key:p/self.pi[key[:3]] for key,p in quartet.items()}
        self.backward={key:p/self.pi[key[1:]] for key,p in quartet.items()}
    def left(self,a,first_two):
        b=a+first_two[0];c=b+first_two[1]
        return [(a-z,self.backward[z,a,b,c]) for z in (a-1,a+1) if 0<=z<self.K]
    def right(self,end,last_two):
        b=end-last_two[-1];a=b-last_two[-2]
        return [(d-end,self.forward[a,b,end,d]) for d in (end-1,end+1) if 0<=d<self.K]
    def cylinder(self,window):
        a,w=window;b=a+w[0]
        if len(w)==1:return self.rho[a,b]
        c=b+w[1];p=self.pi[a,b,c]
        for step in w[2:]:
            d=c+step;p*=self.forward[a,b,c,d];a,b,c=b,c,d
        return p

class FloatMemoryTwo(MemoryTwoLaw):
    def __init__(self,K,up,kappa,eta):
        E=MemoryTwoLaw(K,up,kappa,eta);self.exact=E;self.K=K
        self.rho={key:float(v) for key,v in E.rho.items()}
        self.pi={key:float(v) for key,v in E.pi.items()}
        self.forward={key:float(v) for key,v in E.forward.items()}
        self.backward={key:float(v) for key,v in E.backward.items()}

def moves(model,lev,e):
    r,K,law=model.r,model.K,model.law;assert r>=3
    A,s,B=e[:r],e[r],e[r+1:];end=lev+sum(e)
    left=[(f2,f1,p1*p2) for f1,p1 in law.left(lev,e[:2])
        for f2,p2 in law.left(lev-f1,(f1,e[0]))]
    right=[(g1,g2,p1*p2) for g1,p1 in law.right(end,e[-2:])
        for g2,p2 in law.right(end+g1,(e[-1],g1))]
    out=[('L',(-1,1),[((lev-f1-f2,(f2,f1)+e),p) for f2,f1,p in left]),
        ('R',(1,-1),[((lev,e+(g1,g2)),p) for g1,g2,p in right])]
    splits=[('S1',(1,-1,1))] if s==1 else [('S2',(1,-1,-1)),('S3',(-1,-1,1))]
    for name,trip in splits:
        target=A+trip+B
        if strip_ok(lev,target,K):out.append((name,(1,1),[((lev,target),Q(1))]))
    if s==1 and (A[-1],B[0])!=(1,1):
        core=A[:-1]+(A[-1]+B[0]+1,)+B[1:]
        assert core[:2]==e[:2] and core[-2:]==e[-2:]
        out.append(('M',(-1,-1),[((lev-f1-f2,(f2,f1)+core+(g1,g2)),p*q)
            for f2,f1,p in left for g1,g2,q in right]))
    return out

def build_memory(K,r,up,kappa,eta,exact_control=False):
    assert r>=3;n=2*r+1
    law=MemoryTwoLaw(K,up,kappa,eta) if exact_control else FloatMemoryTwo(K,up,kappa,eta)
    win=list(windows(K,n));rows=list(windows(K,n+2));wi={key:i for i,(key,_) in enumerate(win)};ri={key:i for i,(key,_) in enumerate(rows)}
    bw=[law.cylinder((key>>n,w)) for key,w in win];br=[law.cylinder((key>>(n+2),w)) for key,w in rows]
    rr=array('i');cc=array('i');vv=array('d');keys=[];vw=[];da=[];db=[];starts=[]
    exactrows=[{} for _ in rows] if exact_control else None
    def add(row,col,val):
        rr.append(row);cc.append(col);vv.append(float(val))
        if exact_control:
            val+=exactrows[row].get(col,Q(0))
            if val:exactrows[row][col]=val
            else:exactrows[row].pop(col,None)
    model=SimpleNamespace(K=K,r=r,law=law)
    for i,(key,w) in enumerate(win):
        starts.append(len(keys))
        for label,d,outs in moves(model,key>>n,w):
            col=len(keys);keys.append((key<<3)|LABELS.index(label));vw.append(i);da.append(d[0]);db.append(d[1])
            for (lev,target),p in outs:
                row=ri[(lev<<(n+2))|code(target)];add(row,col,bw[i]*p/br[row])
    starts.append(len(keys))
    for row,(key,w) in enumerate(rows):
        lev=key>>(n+2);inner=((lev+w[0])<<n)|code(w[1:-1]);i=wi[inner]
        for col in range(starts[i],starts[i+1]):add(row,col,Q(-1) if exact_control else -1.)
    M=csr_matrix((np.asarray(vv),(np.asarray(rr),np.asarray(cc))),shape=(len(rows),len(keys)))
    lookup={key:i for i,key in enumerate(keys)}
    flowperm=np.array([lookup[(reflection_key(key>>3,n,K)<<3)|MIRROR[key&7]] for key in keys],dtype=np.int32)
    rowperm=np.array([ri[reflection_key(key,n+2,K)] for key,_ in rows],dtype=np.int32)
    return SimpleNamespace(K=K,r=r,M=M,keys=np.asarray(keys),vw=np.asarray(vw),
        bw=np.asarray(bw,dtype=float),br=np.asarray(br,dtype=float),da=np.asarray(da),db=np.asarray(db),
        flowperm=flowperm,rowperm=rowperm,exactrows=exactrows,law=law)
