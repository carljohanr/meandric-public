"""Exact lift between Model3 objects without refinement:  (rU0,k0,q0,d0) -> (rU,k,q,d), rU>=rU0, k0<=k<=rU+1, q>=q0, d>=d0."""
import numpy as np
from core import head,Q,json
def rev_bits(a,width):
    a=a.astype(np.uint64,copy=True)
    return head.reverse(a,width)
def lift(mo,xo,mn):
    assert mo.rU<=mn.rU and mo.k<=mn.k<=mn.rU+1 and mo.q<=mn.q and mo.d<=mn.d
    xn=np.zeros(mn.nv)
    rU,rU0,k,k0=mn.rU,mo.rU,mn.k,mo.k
    reps,lookup=head.structure(rU);reps=reps.astype(np.uint64)
    inner=((reps>>np.uint64(rU-rU0))&np.uint64((1<<(2*rU0+1))-1)).astype(np.int64)
    val=xo[mo.lookup[inner]].astype(float)
    xF_old=xo[mo.oF:mo.oF+(1<<k0)]
    if k>k0:
        left=reps>>np.uint64(rU+1);right=rev_bits(reps&np.uint64((1<<rU)-1),rU)
        for side in (left,right):
            for j in range(rU-k+1,rU-k0+1):
                gram=((side>>np.uint64(rU-j-k0))&np.uint64((1<<k0)-1)).astype(np.int64)
                val+=xF_old[gram]
    xn[:mn.nu]=val
    g=np.arange(1<<k)>>(k-k0);xn[mn.oF:mn.oF+(1<<k)]=xF_old[g]
    # E
    if mn.ne:
        seen=set()
        for cut,(j,s) in mn.field.items():
            c=int(cut)
            if j<0 or j in seen:continue
            seen.add(j)
            if mo.q==0:continue
            inn=(c>>(mn.q-mo.q))&((1<<(2*mo.q))-1);oj,os=mo.field[str(inn)]
            xn[mn.oE+j]=0 if oj<0 else xo[mo.oE+oj]*os*s
    # G
    if mn.nd and mo.nd:
        done=set()
        for key,j in mn.emap.items():
            if j in done:continue
            done.add(j);w,tri=divmod(key,8);inn=(w>>(mn.d-mo.d))&((1<<(2*mo.d+1))-1)
            xn[mn.oG+j]=xo[mo.oG+mo.emap[8*inn+tri]]
    return xn
