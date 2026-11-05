"""Numerical prototype: projected stationary transport, not literal cycle sums."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import sys,json,pickle,time,gc,argparse,hashlib,resource
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix,vstack
from scipy.sparse.linalg import LinearOperator,lsmr
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parent.parent
sys.path[:0]=[str(ROOT/'1008_lower_122'),str(ROOT/'971_k14_streaming_preflight')]
from cycle_kernel import geometry,kernel
from stream_builder import build
START=time.monotonic()
def event(stage,**kw):
    print(json.dumps(dict(stage=stage,seconds=time.monotonic()-START,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,**kw)),flush=True)
def restriction(fine_keys,r,coarse_keys,r0):
    n=2*r+1;n0=2*r0+1;e=r-r0;assert e>=0
    key=fine_keys.astype(np.int64);win=key>>3;word=win&((1<<n)-1);lev=win>>n;prefix=word>>(n-e);delta=np.zeros(len(key),dtype=np.int64)
    for bit in range(e):delta+=2*((prefix>>bit)&1)-1
    old=(((lev+delta)<<n0)|((word>>e)&((1<<n0)-1)))*8+(key&7)
    order=np.argsort(coarse_keys);sortedkeys=coarse_keys[order];pos=np.searchsorted(sortedkeys,old);assert np.array_equal(sortedkeys[pos],old)
    return order[pos].astype(np.int32)
def load_fields(c,path):
    z=np.load(path,allow_pickle=False);meta=json.loads(str(z['metadata']));r=int(meta['r']);s=geometry(c,r);beta=s.bw[s.vw];j=s.base/(beta@s.base)
    if meta.get('warm_lift_source'):
        parent=Path(meta['warm_lift_source']);s0,j0,_=load_fields(c,parent);mapping=restriction(s.keys,r,s0.keys,s0.r);j=j0[mapping];del s0,j0,mapping;gc.collect()
    B,d=kernel(s);desc=np.array([({'H':0,'L':1,'R':2}[k],key) for k,key in d],dtype=np.int64);assert np.array_equal(z['descriptors'],desc)
    j+=B@z['coefficients'];assert j.min()>0
    del z,B,d,desc;gc.collect();return s,j,meta

def phi(t):
    v=(1+t)*np.log1p(t)-t;small=abs(t)<1e-4;x=t[small]
    v[small]=x*x*(.5+x*(-1/6+x*(1/12+x*(-1/20+x/30))))
    return v

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fields',type=Path,required=True);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--coarse',default='4,5,6');ap.add_argument('--rounds',type=int,default=3);ap.add_argument('--iterations',type=int,default=2000);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.exists()
    c=pickle.loads(a.source.read_bytes());s,j,meta=load_fields(c,a.fields);r=s.r;beta=s.bw[s.vw];j/=beta@j;records=[];packets=[]
    def entropy(j):
        f=beta*j;nu=np.bincount(s.vw,weights=f,minlength=len(s.bw));logs=np.log(f/nu[s.vw]);mass=f.sum();h=-f@logs/mass
        return float(np.exp(2*h)),h,f,nu,logs,mass
    initial=entropy(j)[0];event('fine_witness',r=r,flows=len(j),mu=initial)
    for cr in map(int,a.coarse.split(',')):
        cs=build(c['K'],cr,c['up'],c['kappa']);mapping=restriction(s.keys,r,cs.keys,cr);bc=cs.bw[cs.vw]
        G=vstack((cs.M,csr_matrix(bc*cs.da),csr_matrix(bc*cs.db),csr_matrix(bc)),format='csr');G2=G.copy();G2.data**=2
        event('coarse_built',r=cr,flows=len(cs.keys),rows=cs.M.shape[0],nnz=cs.M.nnz)
        for family in ('head_only','all_moves'):
            selected=np.flatnonzero((cs.keys&7)<2) if family=='head_only' else np.arange(len(cs.keys))
            A=G[:,selected].tocsr();A2=G2[:,selected].tocsr()
            for turn in range(a.rounds if family=='all_moves' else 1):
                old,h,f,nu,logs,mass=entropy(j);g=np.bincount(mapping,weights=beta*(-logs-h),minlength=len(cs.keys));metric=np.bincount(mapping,weights=beta/j,minlength=len(cs.keys));scale=1/np.sqrt(metric[selected]);gg=g[selected]*scale
                rs=1/np.sqrt(np.maximum(A2@(scale*scale),1e-300))
                def mv(x):return rs*(A@(scale*x))
                def rmv(y):return scale*(A.T@(rs*y))
                D=LinearOperator(A.shape,matvec=mv,rmatvec=rmv,dtype=float)
                sol=lsmr(D.T,gg,atol=1e-13,btol=1e-13,conlim=1e12,maxiter=a.iterations)
                q=gg-D.T@sol[0];direction=np.zeros(len(cs.keys));direction[selected]=scale*q
                direction=(direction+direction[cs.flowperm])/2
                # Refine homogeneous feasibility in density coordinates if needed.
                residual=G@direction
                if np.max(abs(residual))>1e-10:
                    row_scale=1/np.sqrt(np.maximum(A2@np.ones(len(selected)),1e-300))
                    T=A.multiply(row_scale[:,None]).tocsr()
                    repair=lsmr(T,-row_scale*residual,atol=1e-14,btol=1e-14,conlim=1e12,maxiter=a.iterations)
                    direction[selected]+=repair[0];residual=G@direction
                fine=direction[mapping];negative=fine<0;maximum=.98*np.min(-j[negative]/fine[negative]) if negative.any() else 10.
                maximum=min(float(maximum),10.);df=beta*fine;dn=np.bincount(s.vw,weights=df,minlength=len(nu));linear=float(df@(-logs-h))
                def gain(t):
                    rel=t*fine/j;rn=t*dn/nu
                    if np.min(rel)<=-1 or np.min(rn)<=-1:return -np.inf
                    return (t*linear-f@phi(rel)+nu@phi(rn))/(mass+t*df.sum())
                opt=minimize_scalar(lambda x:-gain(x),bounds=(0,maximum),method='bounded',options=dict(xatol=1e-10,maxiter=60))
                step=float(opt.x);increment=gain(step);rr=cs.M@direction;flux=abs(cs.M)@abs(direction)
                rawres=float(np.max(abs(rr)));relative=float(np.max(abs(rr)/np.maximum(flux,1e-100)));drifts=(G@direction)[-3:]
                valid=bool(step*rawres<1e-8 and np.max(abs(step*drifts))<1e-10 and increment>1e-13 and np.min(j+step*fine)>0)
                rec=dict(coarse_r=cr,family=family,round=turn,old_mu=old,proposed_mu=float(np.exp(2*(h+increment))),gain_h=increment,step=step,lsmr_iterations=int(sol[2]),lsmr_stop=int(sol[1]),density_residual=step*rawres,relative_direction_residual=relative,drift_mass_residuals=(step*drifts).tolist(),accepted=valid)
                records.append(rec);event('proposal',**rec)
                if valid:
                    j+=step*fine;packets.append((cr,cs.keys.copy(),step*direction));event('accepted',mu=entropy(j)[0])
                del sol,q,fine,df,dn
                if not valid:break
            del A,A2
        del cs,mapping,G,G2;gc.collect()
    final=entropy(j)[0]
    result=dict(status='NUMERICAL_TRANSPORT_FLOW_PROTOTYPE_NOT_CERTIFIED',fine_r=r,initial_mu=initial,final_mu=final,records=records,seconds=time.monotonic()-START,source_fields=str(a.fields),source_fields_sha256=hashlib.sha256(a.fields.read_bytes()).hexdigest(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),notes='Projected directions use all original coarse stationarity rows plus zero drift and zero mass. Only coarse residuals are checked. No full fine-radius balance or exact replay is performed.')
    arrays=dict(metadata=np.array(json.dumps(result)))
    for i,(cr,keys,d) in enumerate(packets):arrays['keys_'+str(i)]=keys;arrays['direction_'+str(i)]=d;arrays['radius_'+str(i)]=np.array(cr)
    with a.out.open('xb') as stream:np.savez(stream,**arrays)
    a.out.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');event('finished',mu=final,packets=len(packets))
if __name__=='__main__':main()
