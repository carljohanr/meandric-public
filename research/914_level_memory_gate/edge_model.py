"""Exact rational second-order nearest-neighbor bulk path law and literal rows."""
from fractions import Fraction as Q
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'905_winding_band_solver'))
from band_model import Model,strip_ok

class EdgeLaw:
    def __init__(self,K,up,kappa=Q(1)):
        self.K=K;self.up={i:Q(up[i]) for i in range(K)};self.kappa=Q(kappa)
        assert K>=2 and self.kappa>0 and self.up[0]==1 and self.up[K-1]==0
        assert all(0<self.up[i]<1 for i in range(1,K-1))
        self.edges=[(a,b) for a in range(K) for b in(a-1,a+1) if 0<=b<K]
        conduct=[Q(1)]
        for i in range(1,K-1):
            conduct.append(conduct[-1]*self.P(i-1,i,i+1)/self.P(i+1,i,i-1))
        norm=2*sum(conduct)
        self.rho={(a,b):conduct[min(a,b)]/norm for a,b in self.edges}
        for b,c in self.edges:
            assert sum(self.rho[a,b]*self.P(a,b,c) for a in(b-1,b+1) if 0<=a<K)==self.rho[b,c]
            assert sum(self.reverse(b,c,a) for a in(b-1,b+1) if 0<=a<K)==1
        assert sum(self.rho.values())==1
    def P(self,a,b,c):
        assert abs(a-b)==1 and 0<=a<self.K and 0<=b<self.K
        if not 0<=c<self.K or abs(c-b)!=1:return Q(0)
        u=self.up[b]
        if u not in(0,1):
            odds=u/(1-u)*self.kappa**(b-a)
            u=odds/(1+odds)
        return u if c==b+1 else 1-u
    def reverse(self,b,c,a):
        return self.rho[a,b]*self.P(a,b,c)/self.rho[b,c]
    def cylinder(self,window):
        a,w=window;b=a+w[0];p=self.rho[a,b]
        for step in w[1:]:
            c=b+step;p*=self.P(a,b,c);a,b=b,c
        return p
    def left(self,a,first):
        b=a+first
        return [(a-prev,self.reverse(a,b,prev)) for prev in(a-1,a+1) if 0<=prev<self.K]
    def right(self,end,last):
        before=end-last
        return [(c-end,self.P(before,end,c)) for c in(end-1,end+1) if 0<=c<self.K]

class EdgeModel(Model):
    def __init__(self,K,r,up,kappa=Q(1)):
        self.law=EdgeLaw(K,up,kappa)
        super().__init__(K,r,self.law.up)
    def constraints(self):
        raise NotImplementedError("Use edge_system.BandSystem: vertex-law Model.constraints is incompatible with edge memory")
    def moves(self,lev,e):
        r,K=self.r,self.K;A,s,B=e[:r],e[r],e[r+1:];end=lev+sum(e)
        left=[(f2,f1,p1*p2) for f1,p1 in self.law.left(lev,e[0])
              for f2,p2 in self.law.left(lev-f1,f1) if p1*p2]
        right=[(g1,g2,p1*p2) for g1,p1 in self.law.right(end,e[-1])
               for g2,p2 in self.law.right(end+g1,g1) if p1*p2]
        out=[('L',(-1,1),[((lev-f1-f2,(f2,f1)+e),p) for f2,f1,p in left]),
             ('R',(1,-1),[((lev,e+(g1,g2)),p) for g1,g2,p in right])]
        splits=[('S1',(1,-1,1))] if s==1 else [('S2',(1,-1,-1)),('S3',(-1,-1,1))]
        for name,trip in splits:
            target=A+trip+B
            if strip_ok(lev,target,K):out.append((name,(1,1),[((lev,target),Q(1))]))
        if s==1 and (A[-1],B[0])!=(1,1):
            core=A[:-1]+(A[-1]+B[0]+1,)+B[1:]
            out.append(('M',(-1,-1),[((lev-f1-f2,(f2,f1)+core+(g1,g2)),p*q)
                       for f2,f1,p in left for g1,g2,q in right]))
        return out
