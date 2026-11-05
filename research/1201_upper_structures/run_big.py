import sys,json,time,argparse,resource,numpy as np
from core3 import *
from opt import anneal
from cand import save_candidate
ap=argparse.ArgumentParser()
ap.add_argument('--a',type=int);ap.add_argument('--cfg');ap.add_argument('--start');ap.add_argument('--betas');ap.add_argument('--iters',type=int,default=3000)
ap.add_argument('--seconds',type=float,default=2400);ap.add_argument('--out');ap.add_argument('--maxcor',type=int,default=15)
A=ap.parse_args()
rU,k,q,d=eval(A.cfg)
t0=time.monotonic();m=Model3(A.a,rU,k,q,d);print('model',m.nv,m.ng,'setup',round(m.setup_seconds,1),'rss MB',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,flush=True)
x=np.load(A.start);ex,s=m.values(x);print('start mu',s.max()**2,flush=True)
def ckpt(xb,Rb):np.save(A.out+'_ckpt.npy',xb)
xb,mu=anneal(m,x,list(map(float,A.betas.split(','))),iters=A.iters,seconds=A.seconds,log=lambda s:print(s,flush=True),checkpoint=ckpt,maxcor=A.maxcor,ckpt_every=200)
np.save(A.out+'_best.npy',xb)
print('final mu',mu,'seconds',time.monotonic()-t0,'peak MB',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,flush=True)
save_candidate(m,xb,A.out+'_cand.json',extra=dict(numerical_mu=mu,seconds=time.monotonic()-t0))
print('saved',flush=True)
