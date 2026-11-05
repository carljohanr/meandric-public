"""Bounded matrix-free radius enlargement. Numerical candidates, not certificates."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import argparse,ctypes,hashlib,json,resource,signal,time,sys
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
LIB=ctypes.CDLL(str(HERE/('evaluator.dylib' if sys.platform == 'darwin' else 'evaluator.so')))
I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
D=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
LIB.evaluate.argtypes=[ctypes.c_int]*5+[I,I,D,ctypes.c_double,
                      ctypes.c_void_p,D,ctypes.c_void_p]
LIB.evaluate.restype=ctypes.c_double

def structure(r):
    n=1<<(2*r+1);w=np.arange(n,dtype=np.uint32);a=w.copy()
    for shift,mask in ((1,0x55555555),(2,0x33333333),(4,0x0f0f0f0f),(8,0x00ff00ff),(16,0x0000ffff)):
        a=((a>>shift)&mask)|((a&mask)<<shift)
    a>>=32-(2*r+1);rep=np.minimum(w,a);representatives=np.flatnonzero(rep==w).astype(np.int32)
    index=np.zeros(n,dtype=np.int32);index[representatives]=np.arange(len(representatives),dtype=np.int32)
    return representatives,index[rep]

class Evaluator:
    def __init__(self,r,k,q):
        if not (2<=r<=12 and 1<=k<=r and 0<=q<=min(r+1,11)): raise ValueError("invalid model dimensions")
        self.q=q;self.field,self.edges,self.ne=fields(q)
        self.r,self.k=r,k;self.reps,self.lookup=structure(r);self.nu=len(self.reps)
        self.maximum=np.zeros(1);self.calls=0;self.seconds=0.
    def evaluate(self,x,beta,gradient=True,rows=False):
        start=time.monotonic();x=np.ascontiguousarray(x,dtype=np.float64)
        if x.shape!=(self.nu+(1<<self.k)+self.ne,) or not np.all(np.isfinite(x)) or not np.isfinite(beta) or beta<=0 or np.max(np.abs(x))>30: raise ValueError("invalid coordinates or beta")
        g=np.zeros(len(x)) if gradient else None
        raw=np.zeros(1<<(2*self.r+2)) if rows else None
        f=LIB.evaluate(self.r,self.k,self.q,self.nu,self.ne,self.lookup,self.edges,x,beta,
                       None if g is None else g.ctypes.data,self.maximum,
                       None if raw is None else raw.ctypes.data)
        self.calls+=1;self.seconds+=time.monotonic()-start
        if rows:return f,float(self.maximum[0]),g,raw
        return f,float(self.maximum[0]),g

def vertex_lift(candidate,ev):
    r,k=ev.r,ev.k;r0,k0=candidate['r'],candidate['k'];assert r0<=r and k0<=k
    old=np.log(np.asarray(candidate['U'],dtype=float)/max(candidate['U']))
    oldf=np.log([float(Fraction(t)) for t in candidate['F']])
    words=ev.reps.astype(np.uint64)
    values=old[(words>>(r-r0))&((1<<(2*r0+1))-1)]
    if k>k0:
        # Prefix extension of F adds the omitted innermost old grams to U.
        for wordside in (words>>(r+1),reverse(words&((1<<r)-1),r)):
            for j in range(r-k+1,r-k0+1):
                values+=oldf[(wordside>>(r-j-k0))&((1<<k0)-1)]
    return np.r_[values,oldf[np.arange(1<<k)>>(k-k0)]]

def reverse(a,length):
    a=a.astype(np.uint64,copy=True)
    for shift,mask in ((1,0x5555555555555555),(2,0x3333333333333333),(4,0x0f0f0f0f0f0f0f0f),
                       (8,0x00ff00ff00ff00ff),(16,0x0000ffff0000ffff),(32,0x00000000ffffffff)):
        a=((a>>shift)&mask)|((a&mask)<<shift)
    return a>>(64-length)

def fields(q):
    if q==0:return {},np.zeros(1,dtype=np.int32),0
    width=2*q;mask=(1<<width)-1;field={};n=0
    for word in range(mask+1):
        if word in field:continue
        rev=int(f'{word:0{width}b}'[::-1],2)
        orbit=[(word,1),(mask-word,1),(rev,-1),(mask-rev,-1)]
        fixed=any(a==b and s!=t for a,s in orbit for b,t in orbit)
        for a,s in orbit:field[a]=[-1,0] if fixed else [n,s]
        if not fixed:n+=1
    field={str(a):field[a] for a in range(mask+1)}
    edges=np.array([0 if field[str(a)][0]<0 else (field[str(a)][0]+1)*field[str(a)][1] for a in range(mask+1)],dtype=np.int32)
    return field,edges,n

def lift(candidate,ev):
    x=vertex_lift(candidate,ev);q0=candidate.get('q',0)
    if q0>ev.q:raise ValueError('cannot shrink edge cut')
    eta=np.zeros(ev.ne)
    if q0:
        old=np.log([float(Fraction(t)) for t in candidate['E']]);seen=set()
        for cut,(j,s) in ev.field.items():
            if j<0 or j in seen:continue
            seen.add(j);inner=(int(cut)>>(ev.q-q0))&((1<<(2*q0))-1)
            oj,os=candidate['field_lookup'][str(inner)]
            eta[j]=0 if oj<0 else old[oj]*os*s
    x[:ev.nu]-=max(x[:ev.nu])
    return np.r_[x,eta]

class Budget(Exception):pass
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--r',type=int,required=True);ap.add_argument('--k',type=int,default=8)
    ap.add_argument('--q',type=int,required=True);ap.add_argument('--scale-u',type=float,default=1);ap.add_argument('--scale-e',type=float,default=1)
    ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--seconds',type=int,default=480);ap.add_argument('--betas',default='5000,20000,80000,320000')
    ap.add_argument('--iterations',type=int,default=500);ap.add_argument('--maxcor',type=int,default=10)
    args=ap.parse_args();assert 2<=args.r<=12 and 1<=args.k<=args.r and args.seconds<=530
    assert 0<args.seconds<=530 and args.iterations>0 and args.maxcor>0
    assert all(np.isfinite(v) and v>0 for v in (args.scale_u,args.scale_e))
    assert not args.out.exists();start=time.monotonic();signal.alarm(args.seconds+55)
    data=args.source.read_bytes();old=json.loads(data);ev=Evaluator(args.r,args.k,args.q);x=lift(old,ev)
    _,starting,_=ev.evaluate(x,1,False)
    print(json.dumps(dict(event='START',r=args.r,k=args.k,variables=len(x),groups=1<<(2*args.r+2),
                          starting_mu=starting**2,evaluator_seconds=ev.seconds)),flush=True)
    best=[starting,x.copy()];trace=[];iterations=[0];lastprint=[time.monotonic()]
    def checkpoint():
        temporary=args.out.with_suffix('.npz.tmp')
        try:
            with temporary.open('wb') as stream:
                np.savez(stream,x=best[1],r=args.r,k=args.k,q=args.q,best_R=best[0])
            temporary.replace(args.out.with_suffix('.npz'))
        except OSError as error:
            print(json.dumps(dict(event='CHECKPOINT_WRITE_FAILED',error=str(error))),flush=True)
    for beta in map(float,args.betas.split(',')):
        def fun(y):
            if time.monotonic()-start>args.seconds:raise Budget()
            assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<1800*2**20
            f,R,g=ev.evaluate(y,beta)
            if R<best[0]:best[:]=[R,y.copy()]
            if time.monotonic()-lastprint[0]>20:
                print(json.dumps(dict(event='PROGRESS',beta=beta,iterations=iterations[0],calls=ev.calls,
                    mu=R*R,best_mu=best[0]**2,seconds=time.monotonic()-start)),flush=True);lastprint[0]=time.monotonic()
            return f,g
        def callback(y):
            iterations[0]+=1
            if iterations[0]%50==0:checkpoint()
        try:
            scale=np.r_[np.full(ev.nu,args.scale_u),np.ones(1<<ev.k),np.full(ev.ne,args.scale_e)]
            bounds=np.array([(-24,12)]*ev.nu+[(-3,2)]*(1<<ev.k)+[(-3,3)]*ev.ne)/scale[:,None]
            if not np.all((x>=bounds[:,0]*scale)&(x<=bounds[:,1]*scale)):raise ValueError('lift outside bounds')
            def scaled(y):
                f,g=fun(y*scale);return f,g*scale
            sol=minimize(scaled,x/scale,jac=True,method='L-BFGS-B',callback=callback,
                         bounds=bounds,
                         options={'maxiter':args.iterations,'maxcor':args.maxcor,'ftol':1e-14,'gtol':1e-9,'maxls':40})
            x=sol.x*scale
            record=dict(beta=beta,status=str(sol.message),iterations=int(sol.nit),best_mu=best[0]**2,
                        seconds=time.monotonic()-start);trace.append(record);print(json.dumps(record),flush=True)
            # OptimizeResult keeps the L-BFGS work array through hess_inv.
            # Release it before the next phase allocates another large array.
            del sol
        except Budget:
            trace.append(dict(beta=beta,status='bounded time budget',best_mu=best[0]**2));break
        checkpoint()
    x=best[1];F=[Fraction(max(1,round(v*2**40)),2**40) for v in np.exp(x[ev.nu:ev.nu+(1<<ev.k)])]
    orbit=np.maximum(1,np.rint(np.exp(x[:ev.nu]-max(x[:ev.nu]))*2**48)).astype(np.int64)
    U=orbit[ev.lookup].tolist()
    E=[Fraction(max(1,round(v*2**36)),2**36) for v in np.exp(x[ev.nu+(1<<ev.k):])]
    rounded=dict(r=args.r,k=args.k,q=args.q,U=U,F=list(map(str,F)),E=list(map(str,E)),field_lookup=ev.field)
    roundedx=lift(rounded,ev);_,roundedR,_=ev.evaluate(roundedx,1,False)
    result=dict(status='NUMERICAL_CANDIDATE_REQUIRES_INDEPENDENT_EXACT_BOUND',**rounded,
        numerical_R=roundedR,numerical_mu=roundedR**2,starting_mu=starting**2,
        raw_U=len(U),independent_U=ev.nu,rows=1<<(2*args.r+3),trace=trace,
        source=str(args.source),source_sha256=hashlib.sha256(data).hexdigest(),
        evaluator_sha256=hashlib.sha256((HERE/'evaluator.cpp').read_bytes()).hexdigest(),
        search_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scale_u=args.scale_u,scale_e=args.scale_e,betas=args.betas,iteration_limit=args.iterations,
        evaluator_binary_sha256=hashlib.sha256((HERE/('evaluator.dylib' if sys.platform == 'darwin' else 'evaluator.so')).read_bytes()).hexdigest(),
        calls=ev.calls,evaluation_seconds=ev.seconds,seconds=time.monotonic()-start,
        peak_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.out.write_text(json.dumps(result)+'\n');checkpoint()
    print(json.dumps({k:v for k,v in result.items() if k not in ('U','F','E','field_lookup','trace')}),flush=True)
if __name__=='__main__':main()
