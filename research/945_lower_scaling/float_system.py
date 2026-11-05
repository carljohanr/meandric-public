"""Stream float numerical operators from frozen914 move rules; never proof data.
Only CSR arrays and window identities persist; no exact coefficient graph.
"""
import sys
from pathlib import Path
from array import array
from types import SimpleNamespace
import numpy as np
import scipy.sparse as sp
R=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(R/'914_level_memory_gate'),str(R/'905_winding_band_solver')]
from edge_model import EdgeLaw,EdgeModel
from edge_system import LAB,word_code
from band_model import all_windows
class FloatLaw:
 def __init__(self,K,up,kappa):
  self.exact=EdgeLaw(K,up,kappa);E=self.exact
  self.rho={e:float(x) for e,x in E.rho.items()}
  self.p={(a,b,c):float(E.P(a,b,c)) for a,b in E.edges for c in (b-1,b+1) if 0<=c<K}
  self.l={(a,b-a):[(d,float(p)) for d,p in E.left(a,b-a)] for a,b in E.edges}
  self.rt={(b,b-a):[(d,float(p)) for d,p in E.right(b,b-a)] for a,b in E.edges}
 def left(self,a,first):return self.l[a,first]
 def right(self,end,last):return self.rt[end,last]
 def cylinder(self,w):
  a,e=w;b=a+e[0];v=self.rho[a,b]
  for d in e[1:]:c=b+d;v*=self.p[a,b,c];a,b=b,c
  return v
class FloatSystem:
 def __init__(self,K,r,up,kappa):
  self.K,self.r,self.up,self.kappa=K,r,up,kappa;self.law=FloatLaw(K,up,kappa)
  wins=all_windows(2*r+1,K);rows=all_windows(2*r+3,K)
  wi={w:i for i,w in enumerate(wins)};ri={w:i for i,w in enumerate(rows)}
  self.nw,self.nx=len(wins),len(rows);self.bw=np.array([self.law.cylinder(w) for w in wins]);br=np.array([self.law.cylinder(w) for w in rows])
  self.keys=[];vw=[];labs=[];da=[];db=[];starts=[]
  rr=array('i');cc=array('i');vv=array('d')
  for i,(lev,e) in enumerate(wins):
   starts.append(len(vw))
   for name,d,outs in EdgeModel.moves(self,lev,e):
    col=len(vw);vw.append(i);labs.append(LAB.index(name));da.append(d[0]);db.append(d[1]);self.keys.append((lev,word_code(e),name))
    for target,p in outs:
     j=ri[target];rr.append(j);cc.append(col);vv.append(self.bw[i]*float(p)/br[j])
  starts.append(len(vw));self.nv=len(vw)
  for j,(lev,e) in enumerate(rows):
   inner=wi[(lev+e[0],e[1:-1])]
   for col in range(starts[inner],starts[inner+1]):rr.append(j);cc.append(col);vv.append(-1.)
  self.vw=np.asarray(vw,dtype=np.int64);self.lab=np.asarray(labs,dtype=np.int8);self.da=np.asarray(da,dtype=np.int8);self.db=np.asarray(db,dtype=np.int8)
  self.M=sp.csr_matrix((np.frombuffer(vv,dtype=np.float64),(np.frombuffer(rr,dtype=np.int32),np.frombuffer(cc,dtype=np.int32))),shape=(self.nx,self.nv))
  self.model=SimpleNamespace(wins=wins)
 def matrix_float(self):return self.M
 def residual(self,j):
  f=self.bw[self.vw]*j
  return max(float(np.max(abs(self.M@j))),abs(float(self.da@f)),abs(float(self.db@f)))
