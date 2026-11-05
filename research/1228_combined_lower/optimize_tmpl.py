"""Optimize selected template directions on top of a baseline (signed columns allowed).  Numerical only."""
from tmpl import *
from scipy.optimize import minimize,Bounds
import argparse,hashlib,resource
ap=argparse.ArgumentParser();ap.add_argument('--selection',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=300);ap.add_argument('--rounds',type=int,default=3)
ap.add_argument('--warm',type=Path,default=None,help='record npz (family mixed_literal_cycles): reoptimize its directions jointly with the new ones, from its coefficients');ap.add_argument('--top',type=int,default=0,help='use only the best N descriptors of the selection');ap.add_argument('--groups',default='',help='restrict to these template groups');a=ap.parse_args()
assert not a.out.exists();start=time.monotonic();selection=json.loads(a.selection.read_text());source=Path(selection['source']);assert hashlib.file_digest(source.open('rb'),'sha256').hexdigest()==selection['source_sha256']
if a.warm is not None:
    with np.load(a.warm) as z:
        wmeta=json.loads(str(z['metadata']));wdesc=z['descriptors'].copy();wcoef=z['coefficients'].copy()
    wsource=Path(wmeta['source']);wsha=wmeta['source_sha256'];assert hashlib.file_digest(wsource.open('rb'),'sha256').hexdigest()==wsha
    assert hashlib.file_digest(a.warm.open('rb'),'sha256').hexdigest()==selection['source_sha256'],'selection must be screened against the warm record'
    base=SparseFineBaseline(wsource,int(selection['r']))
else:
    base=SparseFineBaseline(source,int(selection['r']))
templates=[Template(d['name'],d['kind'],d['labels'],d['m'],d['h0'],d['dtype'],d['side'],d['split']) for d in selection['templates']]
G=TemplateGenerator(base.K,base.r,base.law,base.parameters,templates=templates)
if a.warm is not None:
    from mixed_literal_cycles import LiteralCycleGenerator as Mixed
    class Combined:
        def __init__(self):self.old=Mixed(base.K,base.r,base.law,base.parameters);self.new=G
        def entries(self,d,exact=False,word=None):
            return (self.old if int(d[0])<BASE_KIND else self.new).entries(d,exact=exact,word=word)
    G=Combined()
gof={t.kind:t.name.split(':')[0] for t in templates}
desc=selection['descriptors'];scores=selection['scores']
if a.groups:
    want=set(a.groups.split(','));pick=[i for i,d in enumerate(desc) if gof[d[0]] in want]
else:pick=list(range(len(desc)))
if a.top:pick=pick[:a.top]       # selection is sorted by score
descriptors=[desc[i] for i in pick]
new_count=len(descriptors)
if a.warm is not None:
    descriptors=[[int(k),int(v)] for k,v in wdesc]+descriptors
problem=build_sparse_problem(base,G,descriptors);B=problem['B'];coo=B.tocoo(copy=False);square=B.multiply(B);incidence=np.diff(B.indptr)
theta=np.zeros(len(descriptors))
if a.warm is not None:theta[:len(wcoef)]=wcoef
best=[None,theta.copy()];calls=0;history=[];optimization_start=time.monotonic()
def emit(stage,**kw):
    rec=dict(stage=stage,seconds=round(time.monotonic()-start,1),peak_mib=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,1));rec.update(kw);print(json.dumps(rec),flush=True)
h0_,_,_=sparse_entropy(base,problem,theta);best[0]=h0_;emit('start_point',mu=float(np.exp(2*h0_)),warm=a.warm is not None)
emit('sparse_deep_problem',r=base.r,baseline_mu=base.initial_mu,changed_flows=len(problem['keys']),changed_windows=len(problem['nu']),directions=B.shape[1],nnz=B.nnz)
# NB: incidence = number of directions touching each changed flow (rows of B); B is rows=flows, cols=directions
rowinc=np.diff(B.tocsr().indptr)
class Stop(Exception):pass
for stage in range(a.rounds):
    current_h,_,delta=sparse_entropy(base,problem,theta);current=problem['j']+delta
    scale=1/np.sqrt(np.maximum(np.asarray(square.T@(problem['beta']/current)).ravel(),1e-300))
    pos=coo.data>0;neg=~pos
    fneg=np.full(B.shape[1],np.inf);fpos=np.full(B.shape[1],np.inf)
    np.minimum.at(fneg,coo.col[pos],current[coo.row[pos]]/(rowinc[coo.row[pos]]*coo.data[pos]))
    np.minimum.at(fpos,coo.col[neg],current[coo.row[neg]]/(rowinc[coo.row[neg]]*(-coo.data[neg])))
    oldtheta=theta.copy();best[:]=[current_h,theta.copy()]
    lo=-.999*fneg/scale;hi=.999*fpos/scale
    def fun(y):
        global calls
        if time.monotonic()-optimization_start>a.seconds:raise Stop()
        trial=oldtheta+y*scale;h,g,_=sparse_entropy(base,problem,trial);calls+=1
        if h>best[0]:best[:]=[h,trial.copy()]
        if calls%50==0:emit('sparse_deep_progress',round=stage,calls=calls,mu=float(np.exp(2*h)))
        return -1e4*(h-current_h),-1e4*g*scale
    try:
        sol=minimize(fun,np.zeros(len(theta)),jac=True,method='L-BFGS-B',bounds=Bounds(lo,hi),options=dict(maxiter=800,ftol=1e-16,gtol=1e-11,maxcor=5,maxls=40));history.append(dict(round=stage,iterations=int(sol.nit),message=str(sol.message)))
    except Stop:
        theta=best[1].copy();history.append(dict(round=stage,message='time limit'));break
    theta=best[1].copy();emit('sparse_deep_round',round=stage,mu=float(np.exp(2*best[0])));del sol
h,g,delta=sparse_entropy(base,problem,theta);df=problem['beta']*delta;label=problem['keys']&7
da=np.array([literal.DRIFT[x][0] for x in LABELS]);db=np.array([literal.DRIFT[x][1] for x in LABELS]);drifts=[float(df@da[label]),float(df@db[label])];assert max(map(abs,drifts))<1e-10
used=sorted({d[0] for d in descriptors if d[0]>=BASE_KIND})
metadata=dict(family='template_trial_cycles',status='NUMERICAL_ADAPTIVE_SPARSE_DEEP_CYCLE_NOT_CERTIFIED',r=base.r,K=base.K,initial_mu=base.initial_mu,mu=float(np.exp(2*h)),entropy=float(h),gain_mu=float(np.exp(2*h))-(float(np.exp(2*h0_)) if a.warm is not None else base.initial_mu),start_mu=float(np.exp(2*h0_)),changed_flows=len(problem['keys']),changed_windows=len(problem['nu']),directions=len(theta),nnz=B.nnz,source=str(wsource if a.warm is not None else source),source_sha256=(wsha if a.warm is not None else selection['source_sha256']),warm=str(a.warm) if a.warm is not None else None,selection_sha256=hashlib.file_digest(a.selection.open('rb'),'sha256').hexdigest(),groups=a.groups or 'all',top=a.top,templates=[d for d in selection['templates'] if d['kind'] in used],minimum_changed_density=float(np.min(problem['j']+delta)),both_drift_increment_residuals=drifts,mass=float(base.mass+df.sum()),calls=calls,history=history,seconds=time.monotonic()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,scope='Untouched fine flows inherit positive baseline; all changed flows checked. Numerical only.')
np.savez(a.out,descriptors=np.asarray(descriptors,dtype=np.int64),coefficients=theta,changed_keys=problem['keys'],density_increment=delta,metadata=np.array(json.dumps(metadata)));a.out.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n');emit('finished',**{k:v for k,v in metadata.items() if k not in('templates',)})
