"""Optimize selected deep cycles using only affected flows and window masses."""
from sparse_fine_cycles import *
from scipy.optimize import minimize,Bounds
import argparse,hashlib,resource
ap=argparse.ArgumentParser();ap.add_argument('--selection',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=180);ap.add_argument('--rounds',type=int,default=3);a=ap.parse_args();assert not a.out.exists();start=time.monotonic();selection=json.loads(a.selection.read_text());source=Path(selection['source']);assert hashlib.file_digest(source.open('rb'),'sha256').hexdigest()==selection['source_sha256'];base=SparseFineBaseline(source,int(selection['r']));G=LiteralCycleGenerator(base.K,base.r,base.law,base.parameters);descriptors=selection['descriptors'];problem=build_sparse_problem(base,G,descriptors);B=problem['B'];coo=B.tocoo(copy=False);square=B.multiply(B);incidence=np.diff(B.indptr);theta=np.zeros(len(descriptors));best=[base.h,theta.copy()];calls=0;history=[];optimization_start=time.monotonic()
def emit(stage,**kw):
 rec=dict(stage=stage,seconds=time.monotonic()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20);rec.update(kw);print(json.dumps(rec),flush=True)
emit('sparse_deep_problem',r=base.r,baseline_mu=base.initial_mu,changed_flows=len(problem['keys']),changed_windows=len(problem['nu']),directions=B.shape[1],nnz=B.nnz)
class Stop(Exception):pass
for stage in range(a.rounds):
 current_h,_,delta=sparse_entropy(base,problem,theta);current=problem['j']+delta;scale=1/np.sqrt(np.maximum(np.asarray(square.T@(problem['beta']/current)).ravel(),1e-300));floor=np.full(B.shape[1],np.inf);np.minimum.at(floor,coo.col,current[coo.row]/(incidence[coo.row]*coo.data));oldtheta=theta.copy();best[:]=[current_h,theta.copy()]
 def fun(y):
  global calls
  if time.monotonic()-optimization_start>a.seconds:raise Stop()
  trial=oldtheta+y*scale;h,g,_=sparse_entropy(base,problem,trial);calls+=1
  if h>best[0]:best[:]=[h,trial.copy()]
  if calls%50==0:emit('sparse_deep_progress',round=stage,calls=calls,mu=float(np.exp(2*h)))
  return -1e4*(h-current_h),-1e4*g*scale
 try:
  sol=minimize(fun,np.zeros(len(theta)),jac=True,method='L-BFGS-B',bounds=Bounds(-.999*floor/scale,np.inf),options=dict(maxiter=800,ftol=1e-16,gtol=1e-11,maxcor=5,maxls=40));history.append(dict(round=stage,iterations=int(sol.nit),message=str(sol.message)))
 except Stop:theta=best[1].copy();break
 theta=best[1].copy();emit('sparse_deep_round',round=stage,mu=float(np.exp(2*best[0])));del sol
h,g,delta=sparse_entropy(base,problem,theta);df=problem['beta']*delta;label=problem['keys']&7;da=np.array([literal.DRIFT[x][0] for x in LABELS]);db=np.array([literal.DRIFT[x][1] for x in LABELS]);drifts=[float(df@da[label]),float(df@db[label])];assert max(map(abs,drifts))<1e-10
metadata=dict(status='NUMERICAL_ADAPTIVE_SPARSE_DEEP_CYCLE_NOT_CERTIFIED',r=base.r,K=base.K,initial_mu=base.initial_mu,mu=float(np.exp(2*h)),entropy=float(h),changed_flows=len(problem['keys']),changed_windows=len(problem['nu']),directions=len(theta),nnz=B.nnz,source=str(source),source_sha256=selection['source_sha256'],selection_sha256=hashlib.file_digest(a.selection.open('rb'),'sha256').hexdigest(),minimum_changed_density=float(np.min(problem['j']+delta)),both_drift_increment_residuals=drifts,mass=float(base.mass+df.sum()),calls=calls,history=history,seconds=time.monotonic()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,scope='Untouched fine flows inherit positive baseline; all changed flows checked. Literal cycles passed small exact original-row/drift controls. No full deep stationarity matrix or full deep replay.')
np.savez(a.out,descriptors=np.asarray(descriptors,dtype=np.int64),coefficients=theta,changed_keys=problem['keys'],density_increment=delta,metadata=np.array(json.dumps(metadata)));a.out.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n');emit('finished',**metadata)
