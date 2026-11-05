#!/usr/bin/env python3
"""Independent literal exact verifier for 914-strip-edge-F-v1 pickle certificates.

Reads only the supplied cert. It shares no transition, matrix, or solver code
with 805/905. Mathematical statement and proof are in PROOF.md.
"""
import itertools, pickle, sys, time
from math import lcm
from array import array
import numpy as np
from fractions import Fraction as Q
from pathlib import Path
from mpmath import iv, mp

LABELS=('L','R','S1','S2','S3','M')
DRIFT={'L':(-1,1),'R':(1,-1),'S1':(1,1),'S2':(1,1),'S3':(1,1),'M':(-1,-1)}

def code(word):
    n=0
    for x in word: n=(n<<1)|(x==1)
    return n

def decode(n,length):
    assert isinstance(n,int) and 0<=n<1<<length
    return tuple(1 if (n>>(length-1-i))&1 else -1 for i in range(length))

def admissible(start,word,K):
    a=start
    if not 0<=a<K: return False
    for x in word:
        a+=x
        if not 0<=a<K: return False
    return True

def windows(K,n):
    for a in range(K):
        for w in itertools.product((1,-1),repeat=n):
            if admissible(a,w,K): yield a,w

def literal_move(word,h,label):
    assert 0<h<len(word)-1
    if label=='L':return word,h-1
    if label=='R':return word,h+1
    s=word[h]
    if label=='S1' and s==1:trip=(1,-1,1)
    elif label=='S2' and s==-1:trip=(1,-1,-1)
    elif label=='S3' and s==-1:trip=(-1,-1,1)
    elif label=='M' and s==1 and (word[h-1],word[h+1])!=(1,1):
        x=word[h-1]+word[h+1]+1
        return word[:h-1]+(x,)+word[h+2:],h-1
    else:raise ValueError(f'inapplicable {label} at {word}, head {h}')
    return word[:h]+trip+word[h+1:],h+1

def legal_labels(a,w,r,K):
    out=[]
    for label in LABELS:
        try: new,_=literal_move(w,r,label)
        except ValueError: continue
        if admissible(a,new,K):out.append(label)
    return out

def endpoint(q):
    sign,man,e,_=q._mpf_; n=(-1 if sign else 1)*man
    return Q(n<<e) if e>=0 else Q(n,1<<(-e))

def lower(x):return endpoint(mp.make_mpf(x._mpi_[0]))
def upper(x):return endpoint(mp.make_mpf(x._mpi_[1]))
def I(x):return iv.mpf(x.numerator)/iv.mpf(x.denominator)
def decimal_floor(x,d=12):return Q((x*10**d).__floor__(),10**d)
def decimal_ceil(x,d=12):return Q(-((-x*10**d).__floor__()),10**d)

def parse(cert):
    assert cert['proof_version']=='914-strip-edge-F-v1'
    H,K,r,anchor=(cert[z] for z in ('H','K','r','root_start'))
    assert all(type(v)==int for v in (H,K,r,anchor))
    assert H>=1 and K==2*H and r>=2 and anchor==H-1
    assert len(cert['up'])==K
    up=[]
    for pair in cert['up']:
        assert len(pair)==2 and all(type(t)==int for t in pair) and pair[1]>0
        up.append(Q(*pair))
    assert up[0]==1 and up[-1]==0 and all(0<u<1 for u in up[1:-1])
    den=cert['den'];flows=cert['flows']
    assert type(den)==int and den>0 and isinstance(flows,dict)
    kappa=Q(*cert['kappa']);assert kappa>0
    return H,K,r,up,kappa,den,flows

def stationary_law(up,kappa):
    """Independent exact stationary solve on directed edges (no conductance formula)."""
    K=len(up)
    edges=[(a,b) for a in range(K) for b in range(K) if abs(a-b)==1]
    def P(a,b,c):
        if not(0<=c<K and abs(c-b)==1):return Q(0)
        x=up[b]*kappa**(b-a);u=x/(1-up[b]+x)
        return u if c>b else 1-u
    N=len(edges);rows=[]
    for c,d in edges:
        rows.append([Q(a==c and b==d)-Q(b==c)*P(a,b,d) for a,b in edges]+[Q(0)])
    rows.append([Q(1)]*N+[Q(1)])
    pivotrow=0
    for col in range(N):
        j=next(j for j in range(pivotrow,len(rows)) if rows[j][col])
        rows[pivotrow],rows[j]=rows[j],rows[pivotrow]
        v=rows[pivotrow][col];rows[pivotrow]=[x/v for x in rows[pivotrow]]
        for j in range(len(rows)):
            if j!=pivotrow and rows[j][col]:
                f=rows[j][col];rows[j]=[x-f*y for x,y in zip(rows[j],rows[pivotrow])]
        pivotrow+=1
    rho={edge:rows[i][-1] for i,edge in enumerate(edges)}
    assert sum(rho.values())==1 and all(x>0 for x in rho.values())
    for c,d in edges:
        assert sum(rho[a,b]*P(a,b,d) for a,b in edges if b==c)==rho[c,d]
    return rho,P

def verify(path):
    t0=time.monotonic()
    with Path(path).open('rb') as f:cert=pickle.load(f)
    H,K,r,up,kappa,den,flows=parse(cert)
    rho,P=stationary_law(up,kappa);n=2*r+1
    ws=list(windows(K,n)); rs=list(windows(K,n+2))
    expected={}; nu={}; total=0;drift=[0,0];parity=[0,0]
    for a,w in ws:
        c=code(w);labels=legal_labels(a,w,r,K);s=0
        for m in labels:
            k=(a,c,m); expected[k]=True
            q=flows[k];assert type(q)==int and q>=0,(k,q)
            if m=='L':assert q>0,k
            s+=q;da,db=DRIFT[m];drift[0]+=da*q;drift[1]+=db*q
        nu[a,c]=s;total+=s
        parity[(a+sum(w[:r]))&1]+=s
    assert set(flows)==set(expected),('key mismatch',len(flows),len(expected))
    assert total==den and drift==[0,0] and parity==[den//2,den//2] and den%2==0
    print(f'V1/V3: K={K} r={r}, {len(ws)} windows, {len(flows)} legal flows, '
          f'{len(rs)} original rows; total/drifts/head parity exact',flush=True)

    def right(prev,i):return [(j,j-i,P(prev,i,j)) for j in(i-1,i+1) if 0<=j<K and P(prev,i,j)>0]
    def left(i,nxt):return [(j,i-j,rho[j,i]*P(j,i,nxt)/rho[i,nxt]) for j in(i-1,i+1) if 0<=j<K]
    # Every exterior transition is a rational from a finite directed-edge law.
    # Choose one denominator L for these transitions; four-step products then
    # have denominator dividing D=L**4. Store scaled exact integers throughout.
    right_cache={};left_cache={};L=1
    for i,nxt in rho:
        left_cache[i,nxt]=left(i,nxt)
        right_cache[i,nxt]=right(i,nxt)
        for _,_,p in left_cache[i,nxt]+right_cache[i,nxt]:
            L=lcm(L,p.denominator)
    D=L**4
    left_pairs={};right_pairs={}
    for i,nxt in rho:
        left_pairs[i,nxt]=[(j2,f2,f1,int(p1*p2*L**2))
            for j1,f1,p1 in left_cache[i,nxt]
            for j2,f2,p2 in left_cache[j1,i]]
        right_pairs[i,nxt]=[(k2,g1,g2,int(q1*q2*L**2))
            for k1,g1,q1 in right_cache[i,nxt]
            for k2,g2,q2 in right_cache[nxt,k1]]
        # Prove every conversion above is exact, never truncating a rational.
        assert all((p1*p2*L**2).denominator==1
            for j1,_,p1 in left_cache[i,nxt]
            for _,_,p2 in left_cache[j1,i])
        assert all((q1*q2*L**2).denominator==1
            for k1,_,q1 in right_cache[i,nxt]
            for _,_,q2 in right_cache[nxt,k1])
    print(f'integer scaling: L_bits={L.bit_length()} D_bits={D.bit_length()}',flush=True)
    # Store the small literal incidence graph; accumulate one large-integer
    # row at a time instead of retaining two full vectors of giant integers.
    row_keys=[(a,code(w)) for a,w in rs]
    row_index={key:i for i,key in enumerate(row_keys)}
    flow_values=list(flows.values());flow_index={key:i for i,key in enumerate(flows)}
    graph_rows=array('I');graph_flows=array('I');graph_coeffs=array('I')
    coefficients=[];coefficient_index={}
    for a,w in ws:
        c=code(w);end=a+sum(w)
        for m in legal_labels(a,w,r,K):
            if not flows[a,c,m]:continue
            fi=flow_index[a,c,m]
            for j2,f2,f1,lp in left_pairs[a,a+w[0]]:
                for k2,g1,g2,rp in right_pairs[end-w[-1],end]:
                    V=(f2,f1)+w+(g1,g2);assert admissible(j2,V,K)
                    new,head=literal_move(V,r+2,m)
                    assert admissible(j2,new,K)
                    low=head-r-1;U=new[low:head+r+2]
                    assert len(U)==n+2
                    level=j2+sum(new[:low]);key=(level,code(U))
                    coefficient=lp*rp
                    ci=coefficient_index.get(coefficient)
                    if ci is None:
                        ci=len(coefficients);coefficients.append(coefficient)
                        coefficient_index[coefficient]=ci
                    graph_rows.append(row_index[key]);graph_flows.append(fi);graph_coeffs.append(ci)
    order=np.argsort(np.frombuffer(graph_rows,dtype=np.uint32),kind='stable')
    counts=np.bincount(np.frombuffer(graph_rows,dtype=np.uint32),minlength=len(rs))
    sorted_flows=np.frombuffer(graph_flows,dtype=np.uint32)[order]
    sorted_coeffs=np.frombuffer(graph_coeffs,dtype=np.uint32)[order]
    # memoryview returns Python integer indices without numpy scalar overhead.
    gf=memoryview(sorted_flows);gc=memoryview(sorted_coeffs)
    graph_terms=len(graph_rows)
    del graph_rows,graph_flows,graph_coeffs,order,flow_index,row_index
    print(f'literal graph: terms={graph_terms} coefficients={len(coefficients)}',flush=True)
    offset=0;lhs_total=rhs_total=0
    for ri,(a,w) in enumerate(rs):
        stop=offset+int(counts[ri])
        lhs=sum(flow_values[gf[t]]*coefficients[gc[t]] for t in range(offset,stop))
        inner=(a+w[0],code(w[1:-1]));endinner=inner[0]+sum(w[1:-1])
        coeff=(rho[a,inner[0]]*P(a,inner[0],inner[0]+w[1])/rho[inner[0],inner[0]+w[1]])*P(endinner-w[-2],endinner,endinner+w[-1])*D
        assert coeff.denominator==1
        rhs=nu[inner]*coeff.numerator
        assert lhs==rhs,('stationarity row',row_keys[ri],lhs,rhs)
        lhs_total+=lhs;rhs_total+=rhs;offset=stop
    assert offset==graph_terms
    assert lhs_total==den*D and rhs_total==den*D
    print(f'V2: all {len(rs)} literal original stationarity rows exact',flush=True)
    del sorted_flows,sorted_coeffs,gf,gc,flow_values,counts,row_keys

    iv.prec=mp.prec=200
    ent=iv.mpf(0)
    for (a,c,m),val in flows.items():
        if val:
            ent+=(iv.mpf(val)/iv.mpf(den))*iv.log(iv.mpf(nu[a,c])/iv.mpf(val))
    h_lo=lower(ent);mu_lo=lower(iv.exp(2*ent));mu_hi=upper(iv.exp(2*ent))
    if 'entropy_floor' in cert:
        floor=Q(*cert['entropy_floor']);assert h_lo>=floor,('entropy floor',h_lo,floor)
    if 'mu_floor' in cert:
        floor=Q(*cert['mu_floor']);assert mu_lo>=floor,('mu floor',mu_lo,floor)
    result=dict(H=H,K=K,r=r,root_start=H-1,kappa=str(kappa),windows=len(ws),flows=len(flows),rows=len(rs),
                entropy_lower_12=str(decimal_floor(h_lo)),
                band_mu_lower_12=str(decimal_floor(mu_lo)),
                band_mu_upper_12=str(decimal_ceil(mu_hi)),
                elapsed_seconds=round(time.monotonic()-t0,3))
    print('V4:',result,flush=True)
    return result

if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('usage: verify.py cert.pkl')
    verify(sys.argv[1])
