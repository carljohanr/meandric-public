"""Factor the original head density block as two integer incidences."""
from pathlib import Path
import sys
from collections import deque
from fractions import Fraction as Q
import numpy as np
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'971_k14_streaming_preflight'))
from stream_builder import windows,code

def factor(s):
    K,r=s.K,s.r;n=2*r+1
    heads=np.flatnonzero((s.keys&7)<2)
    headmap={int(s.keys[col]):j for j,col in enumerate(heads)}
    nodes=[key for key,w in windows(K,n+1)];ni={key:i for i,key in enumerate(nodes)}
    rr=[];cc=[];vv=[];component=[]
    for row,key in enumerate(nodes):
        lev=key>>(n+1);word=key&((1<<(n+1))-1);first=1 if word>>n&1 else -1
        left=((lev<<n)|(word>>1))*8+1
        right=(((lev+first)<<n)|(word&((1<<n)-1)))*8
        rr.extend((row,row));cc.extend((headmap[right],headmap[left]));vv.extend((1,-1))
    A=csr_matrix((np.asarray(vv,dtype=np.int64),(rr,cc)),shape=(len(nodes),len(heads)))
    prefix=[];suffix=[]
    for key,w in windows(K,n+2):
        lev=key>>(n+2);word=key&((1<<(n+2))-1);first=w[0]
        prefix.append(ni[(lev<<(n+1))|(word>>1)])
        suffix.append(ni[((lev+first)<<(n+1))|(word&((1<<(n+1))-1))])
    prefix=np.asarray(prefix);suffix=np.asarray(suffix);row=np.repeat(np.arange(len(prefix)),2)
    D=csr_matrix((np.tile(np.array([-1,1],dtype=np.int64),len(prefix)),
        (row,np.column_stack((prefix,suffix)).ravel())),shape=(len(prefix),len(nodes)))
    return heads,A,D,prefix,suffix

def tree_potential(prefix,suffix,forcing):
    """Solve graph gradients by exact tree integration; report all cycle defects."""
    N=int(max(prefix.max(),suffix.max()))+1;adj=[[] for _ in range(N)]
    for row,(u,v) in enumerate(zip(prefix,suffix)):
        adj[u].append((v,row,1));adj[v].append((u,row,-1))
    values=[None]*N;components=0
    for root in range(N):
        if values[root] is not None:continue
        components+=1;values[root]=Q(0);todo=deque([root])
        while todo:
            u=todo.popleft()
            for v,row,sign in adj[u]:
                if values[v] is None:values[v]=values[u]+sign*forcing[row];todo.append(v)
    defects=[values[v]-values[u]-forcing[row] for row,(u,v) in enumerate(zip(prefix,suffix))]
    return values,defects,components

def recover_heads(A,potential):
    """Integrate each small bipartite cut component; return exact defects."""
    N=A.shape[1];adj=[[] for _ in range(N)]
    for row in range(A.shape[0]):
        lo,hi=A.indptr[row:row+2];cols=A.indices[lo:hi];vals=A.data[lo:hi]
        assert len(cols)==2
        plus=int(cols[np.flatnonzero(vals==1)[0]]);minus=int(cols[np.flatnonzero(vals==-1)[0]])
        adj[minus].append((plus,row,1));adj[plus].append((minus,row,-1))
    values=[None]*N;components=0
    for root in range(N):
        if values[root] is not None:continue
        components+=1;values[root]=Q(0);todo=deque([root])
        while todo:
            u=todo.popleft()
            for v,row,sign in adj[u]:
                if values[v] is None:values[v]=values[u]+sign*potential[row];todo.append(v)
    defects=[]
    for row in range(A.shape[0]):
        lo,hi=A.indptr[row:row+2]
        defects.append(sum(int(v)*values[int(col)] for col,v in zip(A.indices[lo:hi],A.data[lo:hi]))-potential[row])
    return values,defects,components
