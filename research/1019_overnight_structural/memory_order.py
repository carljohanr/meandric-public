"""Finite-order reversible path laws with successively preserved marginals."""
from memory_two import *
from memory_two import MemoryTwoLaw

def paths(K,length):
    for a in range(K):
        todo=[(a,)]
        for _ in range(length):todo=[(*v,z) for v in todo for z in (v[-1]-1,v[-1]+1) if 0<=z<K]
        yield from todo

def window(vertices):return vertices[0],tuple(b-a for a,b in zip(vertices,vertices[1:]))

class ExtensionLaw:
    def __init__(self,base,memory,parameter):
        self.base=base;self.K=K=base.K;self.memory=memory;parameter=Q(parameter);assert abs(parameter)<1
        self.state={v:base.cylinder(window(v)) for v in paths(K,memory)}
        joint={v:base.cylinder(window(v)) for v in paths(K,memory+1)}
        def group(v):
            reflected=tuple(K-1-x for x in v);return min(v,v[::-1],reflected,reflected[::-1])
        minimum={}
        for v,p in joint.items():
            if 0<v[1]<K-1 and 0<v[-2]<K-1:
                g=group(v[1:-1]);minimum[g]=min(minimum.get(g,p),p)
        for v in joint:
            if 0<v[1]<K-1 and 0<v[-2]<K-1:joint[v]+=parameter*minimum[group(v[1:-1])]*(1 if v[0]<v[1] else -1)*(1 if v[-1]<v[-2] else -1)
        for v,p in joint.items():assert p>0 and joint[v[::-1]]==p and joint[tuple(K-1-x for x in v)]==p
        for v,p in self.state.items():
            assert sum(joint[(*v,z)] for z in (v[-1]-1,v[-1]+1) if 0<=z<K)==p
            assert sum(joint[(z,*v)] for z in (v[0]-1,v[0]+1) if 0<=z<K)==p
        self.joint=joint;self.forward={v:p/self.state[v[:-1]] for v,p in joint.items()};self.backward={v:p/self.state[v[1:]] for v,p in joint.items()}
    def cylinder(self,w):
        a,steps=w
        if len(steps)<self.memory:return self.base.cylinder(w)
        vertices=[a]
        for step in steps:vertices.append(vertices[-1]+step)
        p=self.state[tuple(vertices[:self.memory+1])]
        for i in range(self.memory+1,len(vertices)):p*=self.forward[tuple(vertices[i-self.memory-1:i+1])]
        return p
    def left(self,a,first):
        assert len(first)==self.memory;vertices=[a]
        for step in first:vertices.append(vertices[-1]+step)
        return [(a-z,self.backward[(z,*vertices)]) for z in (a-1,a+1) if 0<=z<self.K]
    def right(self,end,last):
        assert len(last)==self.memory;vertices=[end]
        for step in last[::-1]:vertices.append(vertices[-1]-step)
        vertices.reverse();return [(z-end,self.forward[(*vertices,z)]) for z in (end-1,end+1) if 0<=z<self.K]

def exact_law(K,up,kappa,parameters):
    assert len(parameters)>=2;law=MemoryTwoLaw(K,up,kappa,parameters[0])
    for memory,parameter in enumerate(parameters[1:],3):law=ExtensionLaw(law,memory,parameter)
    return law

class FloatOrderLaw(ExtensionLaw):
    def __init__(self,K,up,kappa,parameters):
        E=exact_law(K,up,kappa,parameters);self.exact=E;self.K=K;self.memory=E.memory;self.base=E.base
        for name in ('state','forward','backward'):setattr(self,name,{v:float(p) for v,p in getattr(E,name).items()})

def moves(model,lev,e):
    r,K,law=model.r,model.K,model.law;m=law.memory;assert r>=m+1
    A,s,B=e[:r],e[r],e[r+1:];end=lev+sum(e)
    left=[(f2,f1,p1*p2) for f1,p1 in law.left(lev,e[:m])
        for f2,p2 in law.left(lev-f1,(f1,)+e[:m-1])]
    right=[(g1,g2,p1*p2) for g1,p1 in law.right(end,e[-m:])
        for g2,p2 in law.right(end+g1,e[-m+1:]+(g1,))]
    out=[('L',(-1,1),[((lev-f1-f2,(f2,f1)+e),p) for f2,f1,p in left]),
        ('R',(1,-1),[((lev,e+(g1,g2)),p) for g1,g2,p in right])]
    splits=[('S1',(1,-1,1))] if s==1 else [('S2',(1,-1,-1)),('S3',(-1,-1,1))]
    for name,trip in splits:
        target=A+trip+B
        if strip_ok(lev,target,K):out.append((name,(1,1),[((lev,target),Q(1))]))
    if s==1 and (A[-1],B[0])!=(1,1):
        core=A[:-1]+(A[-1]+B[0]+1,)+B[1:]
        assert core[:m]==e[:m] and core[-m:]==e[-m:]
        out.append(('M',(-1,-1),[((lev-f1-f2,(f2,f1)+core+(g1,g2)),p*q)
            for f2,f1,p in left for g1,g2,q in right]))
    return out

def build_memory(K,r,up,kappa,parameters,exact_control=False):
    assert r>=len(parameters)+2;n=2*r+1
    law=exact_law(K,up,kappa,parameters) if exact_control else FloatOrderLaw(K,up,kappa,parameters)
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
