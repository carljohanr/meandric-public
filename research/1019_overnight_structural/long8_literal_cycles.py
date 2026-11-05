"""Sparse positive length-six literal loops with two splits and two merges."""
from sparse_fine_cycles import *
PATTERNS=(('L', 'L', 'S1', 'R', 'R', 'R', 'M', 'L'), ('L', 'L', 'S3', 'R', 'R', 'R', 'M', 'L'), ('L', 'S1', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S1', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S1', 'S3', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S3', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'S3', 'R', 'R', 'M', 'M', 'L'), ('R', 'R', 'S1', 'L', 'L', 'L', 'M', 'R'), ('R', 'R', 'S2', 'L', 'L', 'L', 'M', 'R'), ('R', 'S1', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S1', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S1', 'S3', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S2', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'S3', 'L', 'L', 'M', 'M', 'R'))
class LongCycleGenerator:
 def __init__(self,K,r,law,parameters):
  self.K=K;self.r=r;self.n=2*r+1;self.m=2*r-5;self.law=law;self.parameters=parameters;self.ratios={};assert self.m>=13
 def entries(self,descriptor,exact=False,word=None):
  kind,key=map(int,descriptor);pattern=PATTERNS[kind-3];m=self.m;lev=key>>m;anchor=word if word is not None else literal.decode(key&((1<<m)-1),m);head=m//2;word=anchor;states=[]
  try:
   for label in pattern:
    states.append((word,head,label));word,head=literal.literal_move(word,head,label)
    if not literal.admissible(lev,word,self.K) or word[:3]!=anchor[:3] or word[-3:]!=anchor[-3:]:return []
   assert (word,head)==(anchor,m//2)
  except (ValueError,AssertionError):return []
  begin=m//2-6;local=anchor[begin:begin+13];locallev=lev+sum(anchor[:begin]);cachekey=(exact,locallev,code(local),kind)
  if cachekey not in self.ratios:
   authority=getattr(self.law,'exact',self.law) if exact else self.law;weight=authority.cylinder((locallev,local));ww=local;hh=6;ratios=[]
   for label in pattern:
    ratios.append(weight/authority.cylinder((locallev,ww)));ww,hh=literal.literal_move(ww,hh,label)
   assert (ww,hh)==(local,6);self.ratios[cachekey]=ratios
  out=[]
  for (word,head,label),ratio in zip(states,self.ratios[cachekey]):
   left=self.r-head;right=self.r-(len(word)-1-head);
   if left<0 or right<0:return []
   for exterior in itertools.product((1,-1),repeat=left+right):
    prefix=exterior[:left];full=prefix+word+exterior[left:];start=lev-sum(prefix)
    if literal.admissible(start,full,self.K):out.append(((((start<<self.n)|code(full))<<3)|LABELS.index(label),ratio if exact else float(ratio)))
  return out
LiteralCycleGenerator=LongCycleGenerator
