"""Bounded entropy optimization in the exact local-cycle kernel (numerical stage)."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import argparse,json,pickle,time,resource,hashlib,gc
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,Bounds
from scipy.sparse import csr_matrix
import sys
from local_cycle_family import geometry_witness,kernel_local as kernel
def phi(t):
    v=(1+t)*np.log1p(t)-t;small=abs(t)<1e-4;x=t[small];v[small]=x*x*(.5+x*(-1/6+x*(1/12+x*(-1/20+x/30))));return v

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--r',type=int,required=True);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=300);ap.add_argument('--rounds',type=int,default=8);ap.add_argument('--warm',type=Path);ap.add_argument('--lift-warm',type=Path);ap.add_argument('--objective',choices=['quotient','threshold'],default='quotient');a=ap.parse_args();assert not a.out.exists();start=time.monotonic()
    raw=a.source.read_bytes();s=geometry_witness(a.source,a.r);c=dict(K=s.K);beta=s.bw[s.vw];base=s.base.copy();warm_sha=None
    origin=base.copy()
    B,descriptors=kernel(s);incidence=np.diff(B.indptr);coo=B.tocoo(copy=False);square=csr_matrix((B.data**2,B.indices,B.indptr),shape=B.shape)
    htarget=np.log(12.2)/2;total=np.zeros(B.shape[1]);best=[-np.inf,total.copy(),base.copy()];calls=0
    def score(j):
        f=beta*j;nu=np.bincount(s.vw,weights=f,minlength=len(s.bw));logs=np.log(f/nu[s.vw]);mass=f.sum();H=-f@logs
        value=H/mass if a.objective=='quotient' else H-htarget*mass
        return float(np.exp(2*H/mass)),float(value),logs,mass,H/mass
    if a.warm:
        with np.load(a.warm,allow_pickle=False) as z:
            expected=np.array([({'H':0,'L':1,'R':2}[kind],key) for kind,key in descriptors],dtype=np.int64)
            assert np.array_equal(z['descriptors'],expected);total=z['coefficients'].copy();base+=B@total;assert min(base)>0
    initial=score(base)[0];history=[]
    def event(stage,**kw):
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20
        record=dict(stage=stage,seconds=time.monotonic()-start,peak_mib=peak);record.update(kw)
        print(json.dumps(record),flush=True);assert peak<3800
    event('built',r=a.r,flows=len(base),directions=B.shape[1],nnz=B.nnz,initial_mu=initial)
    setup_seconds=time.monotonic()-start
    optimization_start=time.monotonic()
    class Stop(Exception):pass
    for stage in range(a.rounds):
        floor=np.full(B.shape[1],np.inf);np.minimum.at(floor,coo.col,base[coo.row]/(incidence[coo.row]*coo.data));scale=1/np.sqrt(np.asarray(square.T@(beta/base)).ravel());stagebase=base.copy();stagetotal=total.copy();stagef=beta*stagebase;stagenu=np.bincount(s.vw,weights=stagef,minlength=len(s.bw));stagelogs=np.log(stagef/stagenu[s.vw]);stagemass=stagef.sum();stageh=-stagef@stagelogs/stagemass
        best[:]=[score(base)[1],total.copy(),base.copy()]
        def fun(y):
            nonlocal calls
            if time.monotonic()-optimization_start>a.seconds:raise Stop()
            delta=B@(y*scale);j=stagebase+delta;assert j.min()>0
            df=beta*delta;dn=np.bincount(s.vw,weights=df,minlength=len(s.bw));mass=stagemass+df.sum();stable_gain=(df@(-stagelogs-stageh)-stagef@phi(delta/stagebase)+stagenu@phi(dn/stagenu))/mass;h=stageh+stable_gain;value=h;mu=float(np.exp(2*h));logs=stagelogs+np.log1p(delta/stagebase)-np.log1p(dn/stagenu)[s.vw];g=np.asarray(B.T@(beta*(-logs-(h if a.objective=='quotient' else htarget)))).ravel()*scale
            if a.objective=='quotient':g/=mass
            calls+=1
            if value>best[0]:best[:]=[value,stagetotal+y*scale,j.copy()]
            if calls%50==0:event('progress',round=stage,calls=calls,mu=mu,target_surplus=value)
            # No numerical target cutoff: optimize until budget or convergence.
            return -1e4*stable_gain,-1e4*g
        try:
            sol=minimize(fun,np.zeros(B.shape[1]),jac=True,method='L-BFGS-B',bounds=Bounds(-.999*floor/scale,np.inf),options=dict(maxiter=2000,ftol=1e-16,gtol=1e-12,maxcor=6,maxls=40))
            history.append(dict(round=stage,iterations=int(sol.nit),message=str(sol.message)))
        except Stop:break
        total=best[1].copy();base=best[2].copy();event('round',round=stage,mu=score(base)[0],scaled_gradient_max=float(np.max(np.abs(sol.jac))),near_bound=int(np.count_nonzero(sol.x<=-.999*floor/scale+1e-10)));del sol
    total=best[1];j=best[2];mu,value,*_=score(j)
    # Recompute from the original baseline and saved accumulated coefficients.
    reconstructed=origin+B@total
    assert np.max(abs(reconstructed-j)/j)<1e-10 and reconstructed.min()>0
    metadata=dict(status='NUMERICAL_RAW_MASS_CYCLE_INCREMENT_NOT_CERTIFIED',r=a.r,K=c['K'],mu=mu,initial_mu=initial,target_surplus=value,minimum_density=float(j.min()),source_sha256=hashlib.sha256(raw).hexdigest(),history=history,objective=a.objective,warm_source=str(a.warm) if a.warm else None,warm_lift_source=str(a.lift_warm) if a.lift_warm else None,warm_lift_sha256=warm_sha,setup_seconds=setup_seconds,optimization_seconds=time.monotonic()-optimization_start,calls=calls,seconds=time.monotonic()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20)
    desc=np.array([({'H':0,'L':1,'R':2}[kind],key) for kind,key in descriptors],dtype=np.int64)
    drift_a=float((beta*j)@s.da);drift_b=float((beta*j)@s.db);assert max(abs(drift_a),abs(drift_b))<1e-10
    metadata.update(drift_a=drift_a,drift_b=drift_b,mass=float(beta@j),stationarity_scope='Exact small-cycle comparisons; no full fine-radius stationarity replay. Floating baseline inherited from source.',up=s.source_metadata['up'],kappa=s.source_metadata['kappa'],memory_parameters=list(map(str,s.parameters)),memory_order=len(s.parameters)+1,memory_eta=str(s.parameters[-1]) if s.parameters else None)
    with a.out.open('xb') as f:np.savez(f,keys=s.keys,physical=beta*j,descriptors=desc,coefficients=total,metadata=np.array(json.dumps(metadata)))
    a.out.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n');event('finished',**metadata)
if __name__=='__main__':main()
