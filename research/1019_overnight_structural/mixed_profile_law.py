"""Higher-memory extensions over a height-specific quartet reference law."""
from memory_profile import ProfileLaw,FloatProfileLaw
from memory_order import ExtensionLaw
from fractions import Fraction as Q
class FloatMixedLaw(ExtensionLaw):
 def __init__(self,K,up,kappa,parameters,profile):
  assert len(parameters)>=2;E=ProfileLaw(K,up,kappa,profile)
  for memory,parameter in enumerate(parameters[1:],3):E=ExtensionLaw(E,memory,Q(parameter))
  self.exact=E;self.K=K;self.memory=E.memory;self.base=E.base
  for name in ('state','forward','backward'):setattr(self,name,{key:float(v) for key,v in getattr(E,name).items()})
def make_law(K,up,kappa,parameters,profile):
 return FloatProfileLaw(K,up,kappa,profile) if len(parameters)==1 else FloatMixedLaw(K,up,kappa,parameters,profile)
