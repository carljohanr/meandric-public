"""Sparse exact level-aware strip stationarity in density coordinates.

The enumerator is a frozen copy of 805 flow.py (band_model.py).  Each
equation is the 805 original row divided by its positive row-cylinder mass.
"""
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
from edge_model import EdgeModel

LAB = ['L','R','S1','S2','S3','M']

def word_code(e):
    z=0
    for x in e: z=(z<<1)|(x==1)
    return z

class BandSystem:
    def __init__(self,K,r,up,kappa=F(1)):
        assert K>=2 and r>=2
        self.K,self.r,self.up=K,r,{i:F(up[i]) for i in range(K)}
        assert self.up[0]==1 and self.up[K-1]==0
        assert all(0<self.up[i]<1 for i in range(1,K-1))
        self.kappa=F(kappa)
        self.model=EdgeModel(K,r,self.up,self.kappa)
        self.law=self.model.law
        m=self.model
        self.nw,self.nv,self.nx=len(m.wins),m.nv,len(m.rows)
        self.vw=np.array([wi for wi,_ in m.vars],dtype=np.int64)
        self.lab=np.array([LAB.index(name) for _,name in m.vars],dtype=np.int8)
        self.keys=[(m.wins[wi][0],word_code(m.wins[wi][1]),name) for wi,name in m.vars]
        self.da=np.array([d[0] for d in m.dab],dtype=np.int8)
        self.db=np.array([d[1] for d in m.dab],dtype=np.int8)
        self.be=[self.beta(w) for w in m.wins]
        self.br=[self.beta(w) for w in m.rows]
        self.bw=np.array([float(v) for v in self.be])
        first={}
        for v,(wi,_) in enumerate(m.vars): first.setdefault(wi,[]).append(v)
        rows=[{} for _ in m.rows]
        def add(row,col,val):
            d=rows[row]; v=d.get(col,0)+val
            if v: d[col]=v
            elif col in d: del d[col]
        for v,(wi,_) in enumerate(m.vars):
            for ri,pr in m.trans[v]:
                add(ri,v,self.be[wi]*F(pr)/self.br[ri])
        for ri,(lev,e) in enumerate(m.rows):
            inner=(lev+e[0],e[1:-1]); wi=m.widx[inner]
            for v in first[wi]: add(ri,v,-1)
        self._rows=rows
        # Both drift equations in the original flow measure.
        self._drift=[
            {v:self.be[wi]*int(m.dab[v][j]) for v,(wi,_) in enumerate(m.vars)
             if m.dab[v][j] != 0} for j in (0,1)]
    def beta(self,w): return self.law.cylinder(w)
    def beta_exact(self): return self.be
    def row_dicts(self): return [dict(d) for d in self._rows]+[dict(d) for d in self._drift]
    def matrix_float(self):
        row=[]; col=[]; val=[]
        for i,d in enumerate(self._rows):
            for j,v in d.items(): row.append(i);col.append(j);val.append(float(v))
        return sp.csr_matrix((val,(row,col)),shape=(self.nx,self.nv))
    def residual(self,j):
        A=self.matrix_float(); drift=np.array([sum(float(v)*j[c] for c,v in d.items()) for d in self._drift])
        return max(float(np.max(np.abs(A@j))),float(np.max(np.abs(drift))))
    def shift_only(self):
        j=np.zeros(self.nv)
        j[np.isin(self.lab,[0,1])]=0.5
        return j
    def lift(self,small,js):
        assert self.K==small.K and self.r==small.r+1 and self.up==small.up and self.kappa==small.kappa
        tab={key:v for v,key in enumerate(small.keys)}
        out=np.empty(self.nv)
        for v,(wi,name) in enumerate(self.model.vars):
            lev,e=self.model.wins[wi]; inner=(lev+e[0],e[1:-1])
            key=(inner[0],word_code(inner[1]),name)
            out[v]=js[tab[key]]
        return out
