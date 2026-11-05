import sys,numpy as np
from core3 import *
from cand import save_candidate
a=int(sys.argv[1]);rU,k,q,d=eval(sys.argv[2]);m=Model3(a,rU,k,q,d);x=np.load(sys.argv[3]);ex,s=m.values(x);print('mu',s.max()**2)
save_candidate(m,x,sys.argv[4],extra=dict(numerical_mu=float(s.max()**2)));print('saved')
