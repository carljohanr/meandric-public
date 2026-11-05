"""Independent string/Fraction reference, sharing no producer or C++ implementation."""
from fractions import Fraction as Q

def outgoing(t):
 s,h=t
 if h>0:yield(s,h-1)
 if h+1<len(s):yield(s,h+1)
 for tri in (('101',) if s[h]=='1' else ('100','001')):yield(s[:h]+tri+s[h+1:],h+1)
 if 0<h<len(s)-1:
  rule={'010':'0','011':'1','110':'1'};tri=s[h-1:h+2]
  if tri in rule:yield(s[:h-1]+rule[tri]+s[h+2:],h-1)
def incoming(t):
 s,h=t
 if h>0:yield(s,h-1)
 if h+1<len(s):yield(s,h+1)
 for tri in (('010',) if s[h]=='0' else ('011','110')):yield(s[:h]+tri+s[h+1:],h+1)
 if 0<h<len(s)-1:
  rule={'101':'1','100':'0','001':'0'};tri=s[h-1:h+2]
  if tri in rule:yield(s[:h-1]+rule[tri]+s[h+2:],h-1)
def physical(t):
 s,h=t
 return 0<=h<len(s) and s[0]==s[-1]=='1' and s.count('1')==s.count('0')+1
class Model:
 def __init__(self,data):
  self.r=data['r'];self.k=data['k'];self.q=data['q'];self.F=list(map(Q,data['F']));self.U=data['U'];self.E=list(map(Q,data['E']));self.lookup=data['field_lookup']
 def weight(self,t,adjoint=False,pad=False):
  s,h=t
  if pad:s='0'*(self.r+1)+s+'0'*(self.r+1);h+=self.r+1
  if adjoint:s=s.translate(str.maketrans('01','10'))
  v=Q(self.U[int(s[h-self.r:h+self.r+1],2)])
  for side in (s[:h],s[h+1:][::-1]):
   for j in range(len(side)-self.k+1):v*=self.F[int(side[j:j+self.k],2)]
  return v
 def factor(self,source,target,pad=False):
  s,h=source;v,j=target
  if s!=v or self.q==0:return Q(1)
  assert abs(j-h)==1
  if pad:s='0'*(self.r+1)+s+'0'*(self.r+1);h+=self.r+1
  begin=h-self.q+(j>source[1]);cut=s[begin:begin+2*self.q];assert len(cut)==2*self.q
  index,sign=self.lookup[str(int(cut,2))]
  z=Q(1) if index==-1 else self.E[index]**sign
  return z if j>source[1] else 1/z
 def row(self,t,pad=False,actual=False):
  out=list(outgoing(t));inc=list(incoming(t))
  if actual:out=[a for a in out if physical(a)];inc=[a for a in inc if physical(a)]
  w=self.weight(t,pad=pad);z=self.weight(t,True,pad)
  value=sum((self.weight(a,pad=pad)*self.factor(t,a,pad)/w for a in out),Q(0))
  value+=sum((self.weight(a,True,pad)*self.factor(t,a,pad)/z for a in inc),Q(0))
  return value/2
