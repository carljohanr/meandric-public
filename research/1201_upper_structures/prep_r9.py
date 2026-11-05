import sys,json,time,numpy as np
from core3 import *
from lift3 import lift
from core import head
def x_from_witness3(m,js):
    class EV:pass
    ev=head.Evaluator(m.rU,min(m.k,m.rU),m.q);ev.r=m.rU;ev.k=m.k  # lift uses ev.r, ev.k, reps, ne
    base=head.lift(js,ev)
    x=np.r_[base,np.zeros(m.nd)]
    for key,j in m.emap.items():x[m.oG+j]=np.log(float(Q(js['G'][js['edit_lookup'][str(key)]])))
    return x
if __name__=='__main__':
    src=sys.argv[1];out=sys.argv[2]
    js=json.load(open(src));r,k,q,d=js['r'],js['k'],js['q'],js['d']
    t=time.monotonic();mo=Model3(r+1,r,k,q,d);print('old model',mo.nv,mo.ng,mo.setup_seconds,flush=True)
    xo=x_from_witness3(mo,js);ex,s=mo.values(xo);print('old value',s.max()**2,flush=True)
    r2=r+1;mn=Model3(r2+1,r2,r2+1,r2+1,r2);print('new model',mn.nv,mn.ng,mn.setup_seconds,flush=True)
    xn=lift(mo,xo,mn);ex,s=mn.values(xn);print('new lifted value',s.max()**2,flush=True)
    np.save(out,xn)
