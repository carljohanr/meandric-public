"""Packed original-row numerical builder; no exact graph or BandSystem."""
import sys,time
from pathlib import Path
from array import array
from types import SimpleNamespace
import numpy as np
from scipy.sparse import csr_matrix,vstack
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'945_lower_scaling'))
from float_system import FloatLaw,EdgeModel
LABELS=('L','R','S1','S2','S3','M');MIRROR=(1,0,2,4,3,5)
def decode(c,n):return tuple(1 if c>>(n-1-i)&1 else -1 for i in range(n))
def code(w):
    out=0
    for x in w:out=2*out+(x==1)
    return out
def windows(K,n):
    for a in range(K):
        for c in range((1<<n)-1,-1,-1):
            level=a;w=decode(c,n)
            for x in w:
                level+=x
                if not 0<=level<K:break
            else:yield (a<<n)|c,w
def reflection_key(key,n,K):
    a,c=key>>n,key&((1<<n)-1)
    rev=0
    for i in range(n):rev=2*rev+((c>>i)&1)
    return ((K-1-a-(2*c.bit_count()-n))<<n)|rev
def build(K,r,up,kappa,exact_control=False):
    start=time.monotonic();n=2*r+1;law=FloatLaw(K,up,kappa);exact=law.exact
    assert all(up[K-1-a]==1-up[a] for a in range(K))
    for a,b in exact.edges:
        assert exact.rho[a,b]==exact.rho[b,a]==exact.rho[K-1-a,K-1-b]
        for c in (b-1,b+1):
            if 0<=c<K:
                assert exact.rho[a,b]*exact.P(a,b,c)==exact.rho[c,b]*exact.P(c,b,a)
                assert exact.P(a,b,c)==exact.P(K-1-a,K-1-b,K-1-c)
    if exact_control:law=exact
    win_keys=[];bw=[]
    for key,w in windows(K,n):win_keys.append(key);bw.append(law.cylinder((key>>n,w)))
    wi={key:i for i,key in enumerate(win_keys)}
    row_keys=[];br=[]
    for key,w in windows(K,n+2):row_keys.append(key);br.append(law.cylinder((key>>(n+2),w)))
    ri={key:i for i,key in enumerate(row_keys)};rr=array('i');cc=array('i');vv=array('d')
    exactrows=[{} for _ in row_keys] if exact_control else None
    def add(j,col,value):
        rr.append(j);cc.append(col);vv.append(float(value))
        if exactrows is not None:
            value+=exactrows[j].get(col,0)
            if value:exactrows[j][col]=value
            else:exactrows[j].pop(col,None)
    keys=array('I');vw=array('i');da=array('b');db=array('b');starts=array('i');model=SimpleNamespace(K=K,r=r,law=law)
    for i,key in enumerate(win_keys):
        a=key>>n;w=decode(key&((1<<n)-1),n);starts.append(len(keys))
        for label,d,outs in EdgeModel.moves(model,a,w):
            col=len(keys);keys.append((key<<3)|LABELS.index(label));vw.append(i);da.append(d[0]);db.append(d[1])
            for (lev,target),p in outs:
                j=ri[(lev<<(n+2))|code(target)];add(j,col,bw[i]*(p if exact_control else float(p))/br[j])
    starts.append(len(keys))
    for j,key in enumerate(row_keys):
        a=key>>(n+2);c=key&((1<<(n+2))-1);first=1 if c>>(n+1)&1 else -1
        inner=((a+first)<<n)|((c>>1)&((1<<n)-1));i=wi[inner]
        for col in range(starts[i],starts[i+1]):add(j,col,-1)
    M=csr_matrix((np.frombuffer(vv,dtype=np.float64),(np.frombuffer(rr,dtype=np.int32),np.frombuffer(cc,dtype=np.int32))),shape=(len(row_keys),len(keys)))
    flow_index={key:i for i,key in enumerate(keys)}
    flowperm=np.array([flow_index[(reflection_key(key>>3,n,K)<<3)|MIRROR[key&7]] for key in keys],dtype=np.int32)
    rowperm=np.array([ri[reflection_key(key,n+2,K)] for key in row_keys],dtype=np.int32)
    assert np.array_equal(flowperm[flowperm],np.arange(len(keys))) and np.array_equal(rowperm[rowperm],np.arange(len(row_keys)))
    return SimpleNamespace(M=M,bw=np.asarray(bw,dtype=float),br=np.asarray(br,dtype=float),vw=np.asarray(vw),da=np.asarray(da),db=np.asarray(db),keys=np.asarray(keys),flowperm=flowperm,rowperm=rowperm,row_keys=np.asarray(row_keys),seconds=time.monotonic()-start,K=K,r=r,exactrows=exactrows)
def physical(system):
    beta=system.bw[system.vw]
    return vstack([system.M.multiply(1/beta),csr_matrix(system.da),csr_matrix(system.db)],format='csr'),beta
