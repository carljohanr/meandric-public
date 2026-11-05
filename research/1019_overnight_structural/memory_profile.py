"""Independent exact reversible quartet perturbations by height orbit."""
from memory_two import *
from memory_two import MemoryTwoLaw
class ProfileLaw(MemoryTwoLaw):
 def __init__(self,K,up,kappa,profile):
  base=MemoryTwoLaw(K,up,kappa,Q(0));self.K=K;self.rho=base.rho;self.pi=base.pi;self.edges=base.edges
  quartet=dict(base.quartet);minimum={}
  for (a,b,c,d),p in quartet.items():
   if 0<b<K-1 and 0<c<K-1:
    g=min(b,c,K-1-b,K-1-c);minimum[g]=min(minimum.get(g,p),p)
  self.groups=tuple(sorted(minimum));self.profile={g:Q(profile[g]) for g in self.groups};assert max(map(abs,self.profile.values()))<1
  for key in quartet:
   a,b,c,d=key
   if 0<b<K-1 and 0<c<K-1:
    g=min(b,c,K-1-b,K-1-c);quartet[key]+=self.profile[g]*minimum[g]*(1 if a<b else -1)*(1 if d<c else -1)
  for key,p in quartet.items():assert p>0 and quartet[key[::-1]]==p and quartet[tuple(K-1-z for z in key)]==p
  for (a,b,c),p in self.pi.items():
   assert sum(quartet[a,b,c,d] for d in (c-1,c+1) if 0<=d<K)==p
   assert sum(quartet[z,a,b,c] for z in (a-1,a+1) if 0<=z<K)==p
  self.quartet=quartet;self.forward={key:p/self.pi[key[:3]] for key,p in quartet.items()};self.backward={key:p/self.pi[key[1:]] for key,p in quartet.items()}
class FloatProfileLaw(ProfileLaw):
 def __init__(self,K,up,kappa,profile):
  E=ProfileLaw(K,up,kappa,profile);self.exact=E;self.K=K;self.groups=E.groups;self.profile=E.profile
  for name in ('rho','pi','quartet','forward','backward'):setattr(self,name,{key:float(v) for key,v in getattr(E,name).items()})
