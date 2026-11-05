"""Sparse positive length-six literal loops with two splits and two merges."""
from sparse_fine_cycles import *
PATTERNS=(('L','S1','S3','R','M','M'),('L','S3','S3','R','M','M'),('R','S1','S2','L','M','M'),('R','S2','S2','L','M','M'))
class LongCycleGenerator:
 def __init__(self,K,r,law,parameters):
  self.K=K;self.r=r;self.n=2*r+1;self.m=2*r-5;self.law=law;self.parameters=parameters;self.ratios={};assert self.m>=11
 def entries(self,descriptor,exact=False,word=None):
  kind,key=map(int,descriptor);pattern=PATTERNS[kind-3];m=self.m;lev=key>>m;anchor=word if word is not None else literal.decode(key&((1<<m)-1),m);head=m//2;word=anchor;states=[]
  try:
   for label in pattern:
    states.append((word,head,label));word,head=literal.literal_move(word,head,label)
    if not literal.admissible(lev,word,self.K) or word[:3]!=anchor[:3] or word[-3:]!=anchor[-3:]:return []
   assert (word,head)==(anchor,m//2)
  except (ValueError,AssertionError):return []
  begin=m//2-5;local=anchor[begin:begin+11];locallev=lev+sum(anchor[:begin]);cachekey=(exact,locallev,code(local),kind)
  if cachekey not in self.ratios:
   authority=getattr(self.law,'exact',self.law) if exact else self.law;weight=authority.cylinder((locallev,local));ww=local;hh=5;ratios=[]
   for label in pattern:
    ratios.append(weight/authority.cylinder((locallev,ww)));ww,hh=literal.literal_move(ww,hh,label)
   assert (ww,hh)==(local,5);self.ratios[cachekey]=ratios
  out=[]
  for (word,head,label),ratio in zip(states,self.ratios[cachekey]):
   left=self.r-head;right=self.r-(len(word)-1-head);assert left>=0 and right>=0
   for exterior in itertools.product((1,-1),repeat=left+right):
    prefix=exterior[:left];full=prefix+word+exterior[left:];start=lev-sum(prefix)
    if literal.admissible(start,full,self.K):out.append(((((start<<self.n)|code(full))<<3)|LABELS.index(label),ratio if exact else float(ratio)))
  return out
LiteralCycleGenerator=LongCycleGenerator
