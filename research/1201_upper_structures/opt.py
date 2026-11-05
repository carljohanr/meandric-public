"""Smoothed-max L-BFGS annealing for any model with X, owner, multiplicity (Model/Model2/Model3)."""
import time,json,resource
import numpy as np
from scipy.optimize import minimize
class Budget(Exception):pass
def default_scale(m):
    nf=1<<m.k
    parts=[np.full(m.nu,np.sqrt(m.nu/nf)),np.ones(nf),np.full(m.ne,np.sqrt(max(1,m.ne)/nf)),np.full(m.nd,np.sqrt(max(1,m.nd)/nf))]
    rest=m.nv-sum(len(p) for p in parts)
    if rest>0:parts.append(np.full(rest,np.sqrt(max(1,rest)/nf)))
    return np.concatenate(parts)
def anneal(m,x,betas,iters=1500,seconds=600,log=print,target=0,lo=-24,hi=12,maxcor=15,checkpoint=None,ckpt_every=300):
    ng=m.ng;scale=default_scale(m);start=time.monotonic()
    ex=np.exp(m.X@x);sums=np.bincount(m.owner,weights=ex,minlength=ng)/2
    best=[float(sums.max()),x.copy()];last=[time.monotonic()];nfev=[0]
    lob=np.full(m.nv,-4.0);hib=np.full(m.nv,4.0);lob[:m.nu]=lo;hib[:m.nu]=hi
    # E/G/F bounds like original (-3..3)
    bounds=np.c_[lob,hib]/scale[:,None]
    x=np.clip(x,lob,hib)
    for beta in betas:
        def fun(y):
            if time.monotonic()-start>seconds:raise Budget()
            xx=y*scale;ex=np.exp(m.X@xx);sums=np.bincount(m.owner,weights=ex,minlength=ng)/2;R=float(sums.max());logs=np.log(sums)
            nfev[0]+=1
            if R<best[0]:best[:]=[R,xx.copy()]
            if target and best[0]**2<=target:raise Budget()
            if checkpoint and nfev[0]%ckpt_every==0:checkpoint(best[1],best[0])
            if time.monotonic()-last[0]>30:log(json.dumps(dict(ev='PROG',beta=beta,mu=R*R,best_mu=best[0]**2,t=round(time.monotonic()-start))));last[0]=time.monotonic()
            p=m.multiplicity*np.exp(beta*(logs-logs.max()));tot=p.sum();p/=tot
            g=np.asarray(m.X.T@(p[m.owner]*ex/(2*sums[m.owner]))).ravel()
            return float(logs.max()+np.log(tot)/beta),g*scale
        try:
            sol=minimize(fun,x/scale,jac=True,method='L-BFGS-B',bounds=bounds,options=dict(maxiter=iters,ftol=1e-15,gtol=1e-10,maxls=40,maxcor=maxcor))
            x=sol.x*scale;log(json.dumps(dict(ev='STAGE',beta=beta,nit=int(sol.nit),best_mu=best[0]**2,t=round(time.monotonic()-start))));del sol
        except Budget:break
    if checkpoint:checkpoint(best[1],best[0])
    return best[1],best[0]**2
