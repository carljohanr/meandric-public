"""Matrix-free head elimination; exact graph topology, floating edit coefficients.

The compatibility operator vanishes exactly when edits admit a stationary head
response. Its adjoint allows iterative constrained solves without dense paths.
"""
from collections import deque
import numpy as np
from scipy.sparse.linalg import LinearOperator
from head_factor import factor

class TreeIntegrator:
    def __init__(self,minus,plus,nodes):
        self.minus=np.asarray(minus);self.plus=np.asarray(plus);self.nodes=nodes
        adj=[[] for _ in range(nodes)]
        for edge,(u,v) in enumerate(zip(minus,plus)):
            adj[u].append((v,edge,1));adj[v].append((u,edge,-1))
        parent=np.full(nodes,-1,dtype=int);edge=np.full(nodes,-1,dtype=int);sign=np.zeros(nodes);depth=np.full(nodes,-1,dtype=int);roots=[]
        for root in range(nodes):
            if depth[root]>=0:continue
            roots.append(root);depth[root]=0;todo=deque([root])
            while todo:
                u=todo.popleft()
                for v,e,s in adj[u]:
                    if depth[v]<0:
                        parent[v]=u;edge[v]=e;sign[v]=s;depth[v]=depth[u]+1;todo.append(v)
        self.parent=parent;self.edge=edge;self.sign=sign;self.roots=np.asarray(roots)
        self.levels=[np.flatnonzero(depth==d) for d in range(1,int(depth.max())+1)]
        self.shape=(nodes,len(minus))
    def apply(self,b):
        x=np.zeros(self.nodes)
        for v in self.levels:x[v]=x[self.parent[v]]+self.sign[v]*b[self.edge[v]]
        return x
    def adjoint(self,x):
        acc=np.array(x,copy=True);b=np.zeros(self.shape[1])
        for v in self.levels[::-1]:
            b[self.edge[v]]=self.sign[v]*acc[v]
            np.add.at(acc,self.parent[v],acc[v])
        return b

def compatibility(s):
    heads,A,D,prefix,suffix=factor(s)
    edits=np.flatnonzero((s.keys&7)>=2);E=s.M[:,edits].tocsr()
    td=TreeIntegrator(prefix,suffix,D.shape[1])
    minus=[];plus=[]
    for row in range(A.shape[0]):
        lo,hi=A.indptr[row:row+2];cols=A.indices[lo:hi];vals=A.data[lo:hi]
        minus.append(int(cols[np.flatnonzero(vals==-1)[0]]));plus.append(int(cols[np.flatnonzero(vals==1)[0]]))
    ta=TreeIntegrator(minus,plus,A.shape[1]);m=D.shape[0];n=A.shape[0]
    def matvec(e):
        b=-(E@e);y=td.apply(b);h=ta.apply(y)
        return np.r_[b-D@y,y-A@h]
    def rmatvec(z):
        u,v=z[:m],z[m:]
        yy=v-ta.adjoint(A.T@v)-D.T@u
        return -(E.T@(u+td.adjoint(yy)))
    C=LinearOperator((m+n,len(edits)),matvec=matvec,rmatvec=rmatvec,dtype=float)
    def reconstruct(e):
        return ta.apply(td.apply(-(E@e)))
    return C,heads,edits,reconstruct,td,ta
