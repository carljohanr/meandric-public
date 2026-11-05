from sparse_fine_cycles import *
from sparse_fine_cycles import LiteralCycleGenerator as Short
from long_literal_cycles import LongCycleGenerator as Six
from long8_literal_cycles import LongCycleGenerator as Eight
class LiteralCycleGenerator:
 def __init__(self,K,r,law,parameters):self.short=Short(K,r,law,parameters);self.six=Six(K,r,law,parameters);self.eight=Eight(K,r,law,parameters)
 def entries(self,descriptor,exact=False,word=None):
  kind,key=map(int,descriptor)
  if kind<3:return self.short.entries((kind,key),exact=exact,word=word)
  if kind<7:return self.six.entries((kind,key),exact=exact,word=word)
  return self.eight.entries((kind-4,key),exact=exact,word=word)
