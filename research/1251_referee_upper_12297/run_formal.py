"""Exact formal-row run: all 2^(2r+3) rows (or a slice) of a witness.  usage: run_formal.py cand.json out.json [part nparts]"""
import sys,json,time,hashlib
sys.dont_write_bytecode=True
from fractions import Fraction as Q
import mychk as M
path,outp=sys.argv[1],sys.argv[2]
part,nparts=(int(sys.argv[3]),int(sys.argv[4])) if len(sys.argv)>4 else (0,1)
claimed=sys.argv[5] if len(sys.argv)>5 else None
t0=time.time()
c=M.load(path); r=c.r; L=2*r+3; total=1<<L
Rn,Rd=(Q(json.load(open(claimed))['R']).numerator,Q(json.load(open(claimed))['R']).denominator) if claimed else (None,None)
lo=total*part//nparts; hi=total*(part+1)//nparts
viol=0; best=None; fl=[]
import array
vals=array.array('d')
for code in range(lo,hi):
    N,D=M.row(c,(format(code,'0%db'%L),r+1))
    v=N/(2*D)           # float, correctly rounded big-int division
    vals.append(v)
    if best is None or N*best[1]>best[0]*D: best=(N,D,code)
    if Rn is not None and not M.leq(N,D,Rn,Rd): viol+=1
    if (code-lo)%200000==0: print(code-lo,time.time()-t0,flush=True)
res=dict(r=r,k=c.k,q=c.q,d=c.d,part=part,nparts=nparts,rows=hi-lo,violations_vs_claimed=viol,
   argmax=best[2],Rmax=str(Q(best[0],2*best[1])),Rmax_float=float(Q(best[0],2*best[1])),mu_float=float(Q(best[0],2*best[1])**2),
   seconds=time.time()-t0,cand_sha256=hashlib.sha256(open(path,'rb').read()).hexdigest())
json.dump(res,open(outp,'w'),indent=1)
vals.tofile(open(outp+'.f64','wb'))
print({k:v for k,v in res.items() if k!='Rmax'})
