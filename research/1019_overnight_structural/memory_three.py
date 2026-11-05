"""Positive reversible quintet perturbation preserving all quartet marginals.

Experimental order-three boundary law. Not an accepted growth certificate.
"""
from memory_two import *
from memory_two import MemoryTwoLaw

class MemoryThreeLaw:
    def __init__(self,K,up,kappa,eta,zeta):
        base=MemoryTwoLaw(K,up,kappa,eta);self.K=K;self.eta=Q(eta);self.zeta=Q(zeta);assert abs(self.zeta)<1
        self.rho=base.rho;self.pi=base.pi;self.quartet=base.quartet
        quintet={(*key,e):p*base.forward[(*key[1:],e)] for key,p in self.quartet.items() for e in (key[-1]-1,key[-1]+1) if 0<=e<K}
        def group(middle):
            reverse=middle[::-1];reflected=tuple(K-1-x for x in middle)
            return min(middle,reverse,reflected,reflected[::-1])
        minimum={}
        for key,p in quintet.items():
            a,b,c,d,e=key
            if 0<b<K-1 and 0<d<K-1:
                g=group((b,c,d));minimum[g]=min(minimum.get(g,p),p)
        for key in quintet:
            a,b,c,d,e=key
            if 0<b<K-1 and 0<d<K-1:quintet[key]+=self.zeta*minimum[group((b,c,d))]*(1 if a<b else -1)*(1 if e<d else -1)
        for key,p in quintet.items():
            assert p>0 and quintet[key[::-1]]==p and quintet[tuple(K-1-x for x in key)]==p
        for key,p in self.quartet.items():
            assert sum(quintet[(*key,e)] for e in (key[-1]-1,key[-1]+1) if 0<=e<K)==p
            assert sum(quintet[(a,*key)] for a in (key[0]-1,key[0]+1) if 0<=a<K)==p
        self.quintet=quintet;self.forward={key:p/self.quartet[key[:4]] for key,p in quintet.items()};self.backward={key:p/self.quartet[key[1:]] for key,p in quintet.items()}
    def left(self,a,first_three):
        vertices=[a]
        for step in first_three:vertices.append(vertices[-1]+step)
        return [(a-z,self.backward[(z,*vertices)]) for z in (a-1,a+1) if 0<=z<self.K]
    def right(self,end,last_three):
        vertices=[end]
        for step in last_three[::-1]:vertices.append(vertices[-1]-step)
        vertices.reverse()
        return [(z-end,self.forward[(*vertices,z)]) for z in (end-1,end+1) if 0<=z<self.K]
    def cylinder(self,window):
        a,w=window;vertices=[a]
        for step in w:vertices.append(vertices[-1]+step)
        if len(w)==1:return self.rho[tuple(vertices)]
        if len(w)==2:return self.pi[tuple(vertices)]
        p=self.quartet[tuple(vertices[:4])]
        for i in range(4,len(vertices)):p*=self.forward[tuple(vertices[i-4:i+1])]
        return p

class FloatMemoryThree(MemoryThreeLaw):
    def __init__(self,K,up,kappa,eta,zeta):
        E=MemoryThreeLaw(K,up,kappa,eta,zeta);self.exact=E;self.K=K
        for name in ('rho','pi','quartet','forward','backward'):setattr(self,name,{key:float(v) for key,v in getattr(E,name).items()})

def moves(model,lev,e):
    r,K,law=model.r,model.K,model.law;assert r>=4
    A,s,B=e[:r],e[r],e[r+1:];end=lev+sum(e)
    left=[(f2,f1,p1*p2) for f1,p1 in law.left(lev,e[:3])
        for f2,p2 in law.left(lev-f1,(f1,e[0],e[1]))]
    right=[(g1,g2,p1*p2) for g1,p1 in law.right(end,e[-3:])
        for g2,p2 in law.right(end+g1,(e[-2],e[-1],g1))]
    out=[('L',(-1,1),[((lev-f1-f2,(f2,f1)+e),p) for f2,f1,p in left]),
        ('R',(1,-1),[((lev,e+(g1,g2)),p) for g1,g2,p in right])]
    splits=[('S1',(1,-1,1))] if s==1 else [('S2',(1,-1,-1)),('S3',(-1,-1,1))]
    for name,trip in splits:
        target=A+trip+B
        if strip_ok(lev,target,K):out.append((name,(1,1),[((lev,target),Q(1))]))
    if s==1 and (A[-1],B[0])!=(1,1):
        core=A[:-1]+(A[-1]+B[0]+1,)+B[1:]
        assert core[:3]==e[:3] and core[-3:]==e[-3:]
        out.append(('M',(-1,-1),[((lev-f1-f2,(f2,f1)+core+(g1,g2)),p*q)
            for f2,f1,p in left for g1,g2,q in right]))
    return out

def build_memory(K,r,up,kappa,eta,zeta,exact_control=False):
    assert r>=4;n=2*r+1
    law=MemoryThreeLaw(K,up,kappa,eta,zeta) if exact_control else FloatMemoryThree(K,up,kappa,eta,zeta)
    win=list(windows(K,n));rows=list(windows(K,n+2));wi={key:i for i,(key,_) in enumerate(win)};ri={key:i for i,(key,_) in enumerate(rows)}
    bw=[law.cylinder((key>>n,w)) for key,w in win];br=[law.cylinder((key>>(n+2),w)) for key,w in rows]
    rr=array('i');cc=array('i');vv=array('d');keys=[];vw=[];da=[];db=[];starts=[]
    exactrows=[{} for _ in rows] if exact_control else None
    def add(row,col,val):
        rr.append(row);cc.append(col);vv.append(float(val))
        if exact_control:
            val+=exactrows[row].get(col,Q(0))
            if val:exactrows[row][col]=val
            else:exactrows[row].pop(col,None)
    model=SimpleNamespace(K=K,r=r,law=law)
    for i,(key,w) in enumerate(win):
        starts.append(len(keys))
        for label,d,outs in moves(model,key>>n,w):
            col=len(keys);keys.append((key<<3)|LABELS.index(label));vw.append(i);da.append(d[0]);db.append(d[1])
            for (lev,target),p in outs:
                row=ri[(lev<<(n+2))|code(target)];add(row,col,bw[i]*p/br[row])
    starts.append(len(keys))
    for row,(key,w) in enumerate(rows):
        lev=key>>(n+2);inner=((lev+w[0])<<n)|code(w[1:-1]);i=wi[inner]
        for col in range(starts[i],starts[i+1]):add(row,col,Q(-1) if exact_control else -1.)
    M=csr_matrix((np.asarray(vv),(np.asarray(rr),np.asarray(cc))),shape=(len(rows),len(keys)))
    lookup={key:i for i,key in enumerate(keys)}
    flowperm=np.array([lookup[(reflection_key(key>>3,n,K)<<3)|MIRROR[key&7]] for key in keys],dtype=np.int32)
    rowperm=np.array([ri[reflection_key(key,n+2,K)] for key,_ in rows],dtype=np.int32)
    return SimpleNamespace(K=K,r=r,M=M,keys=np.asarray(keys),vw=np.asarray(vw),
        bw=np.asarray(bw,dtype=float),br=np.asarray(br,dtype=float),da=np.asarray(da),db=np.asarray(db),
        flowperm=flowperm,rowperm=rowperm,exactrows=exactrows,law=law)
