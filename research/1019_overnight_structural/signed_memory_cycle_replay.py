"""Independent cycle-witness verifier: accepted exact baseline + dyadic cycles.

No producer, sparse matrix, optimizer, or numerical candidate density is used.
The local-cycle identity is checked through literal closed paths and unchanged
boundary edges. Every fine legal flow is checked for strict positivity using
directed MPFR enclosures, and conditional entropy is rounded outward.
"""
import sys,json,pickle,hashlib,time,gc,argparse,itertools,importlib.util,resource
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
from memory_three import MemoryThreeLaw
sys.path.insert(0,str(ROOT/'910_band_bound_push/vendor'))
import gmpy2 as g
spec=importlib.util.spec_from_file_location('literal940cycle',ROOT/'940_integer_literal_verifier/verify_sorted.py');literal=importlib.util.module_from_spec(spec);spec.loader.exec_module(literal)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for data in iter(lambda:f.read(8*2**20),b''):h.update(data)
    return h.hexdigest()

def terms(a,c,label,r):
    """Invert each possible position of this edge in a two/four-move cycle."""
    n=2*r+1;m=n-4;ell=n-2;h=r-2
    def level(start):return a+2*(c>>(n-start)).bit_count()-start if start else a
    def short(start):return (level(start)<<m)|((c>>(n-start-m))&((1<<m)-1))
    def long(start,side):
        core=(c>>(n-start-ell))&((1<<ell)-1);trip=(core>>h)&7
        if (side==1 and trip&6!=4) or (side==2 and trip&3!=1):return None
        bit=(trip&1) if side==1 else trip>>2
        anchor=((core>>(h+3))<<(h+1))|(bit<<h)|(core&((1<<h)-1))
        return (level(start)<<m)|anchor
    out=[]
    if label in ('L','R'):
        cut=(a<<(n-1))|(c>>1) if label=='L' else ((a+(1 if c>>(n-1)&1 else -1))<<(n-1))|(c&((1<<(n-1))-1))
        out.append((0,cut,False))
    if label in ('S1','S2'):out.append((1,short(2),False))
    if label in ('S1','S3'):out.append((2,short(2),False))
    if label=='L':
        key=long(1,1)
        if key is not None:out.append((1,key,True))
        out.append((2,short(1),False))
    if label=='R':
        key=long(1,2)
        if key is not None:out.append((2,key,True))
        out.append((1,short(3),False))
    if label=='M':
        for side,start in ((1,2),(2,0)):
            key=long(start,side)
            if key is not None:out.append((side,key,True))
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fields',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--precision',type=int,default=160);ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--parent-replay',type=Path,required=True);a=ap.parse_args();assert not a.out.exists();assert a.precision>=128;start=time.monotonic()
    g.get_context().precision=a.precision
    down=g.context(precision=a.precision,round=g.RoundDown);up=g.context(precision=a.precision,round=g.RoundUp)
    def event(stage,**kw):print(json.dumps(dict(stage=stage,seconds=time.monotonic()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,**kw)),flush=True)
    def enclose(q):
        if isinstance(q,Q):q=g.mpq(q.numerator,q.denominator)
        with down:lo=g.mpfr(q)
        with up:hi=g.mpfr(q)
        return lo,hi
    accepted=a.parent_replay;manifest=json.loads(accepted.read_text());baseline=a.parent
    assert manifest['status'] in ('PASS_INDEPENDENT_ORIGINAL_ROWS_POSITIVITY_DRIFTS_AND_OUTWARD_ENTROPY_PARENT','PASS_INDEPENDENT_ORIGINAL_ROWS_DRIFTS_SIGNED_STATIONARY_PARENT')
    assert sha(baseline)==manifest['certificate_sha256']
    fields_sha=sha(a.fields)
    with np.load(a.fields,allow_pickle=False) as z:descriptors=z['descriptors'];coeff=z['coefficients'];meta=json.loads(str(z['metadata']))
    cert=pickle.load(baseline.open('rb'));K,oldr,lawup,kappa,D=cert['K'],cert['r'],cert['up'],cert['kappa'],cert['denominator'];r=int(meta['r']);n=2*r+1;oldn=2*oldr+1;m=n-4;memory=3
    assert K==meta['K'] and r>=oldr and r>=memory+2 and tuple(map(Q,meta['memory_parameters']))==cert['parameters']
    law=MemoryThreeLaw(K,lawup,kappa,*cert['parameters']);PI={key:enclose(value) for key,value in law.forward.items()}
    def all_windows(width):
        def visit(first,vertices,word,c,lo,hi):
            if len(word)==width:yield first,word,c,lo,hi;return
            current=vertices[-1]
            for nxt in (current+1,current-1):
                if 0<=nxt<K:
                    pl,pu=PI[(*vertices,nxt)]
                    yield from visit(first,(*vertices[1:],nxt),word+(nxt-current,),(c<<1)|(nxt>current),down.mul(lo,pl),up.mul(hi,pu))
        for first,word in literal.windows(K,memory):
            vertices=[first]
            for step in word:vertices.append(vertices[-1]+step)
            lo,hi=enclose(law.cylinder((first,word)))
            yield from visit(first,tuple(vertices),word,literal.code(word),lo,hi)
    def cylinder(lev,w):return law.cylinder((lev,w))
    flows={}
    for packed,value in zip(cert['keys'],cert['numerators']):
        packed=int(packed);window=packed>>3;flows[window>>oldn,window&((1<<oldn)-1),literal.LABELS[packed&7]]=value
    Dlo,Dhi=enclose(D);coarse={};total=0;drift=[0,0];count=0
    for lev,w,c,blo,bhi in all_windows(oldn):
        for label in literal.legal_labels(lev,w,oldr,K):
            N=flows[lev,c,label];assert type(N) is int
            nl,nh=enclose(N)
            if N>=0:lo=down.div(down.div(nl,Dhi),bhi);hi=up.div(up.div(nh,Dlo),blo)
            else:lo=down.div(down.div(nl,Dlo),blo);hi=up.div(up.div(nh,Dhi),bhi)
            coarse[lev,c,label]=(lo,hi);count+=1;total+=N
            dl,dr=literal.DRIFT[label];drift[0]+=dl*N;drift[1]+=dr*N
    assert count==len(flows) and total==D and drift==[0,0]
    del cert,flows;gc.collect();event('baseline',flows=count,accepted_lower=manifest['mu_lower_decimal'])
    fields={};cycle_count=0;head_count=0
    for (kind,key),value in zip(descriptors,coeff):
        kind,key=int(kind),int(key);assert kind in (0,1,2) and np.isfinite(value) and (kind,key) not in fields
        point=Q.from_float(float(value));zl,zh=enclose(point)
        if kind==0:
            word=literal.decode(key&((1<<(n-1))-1),n-1);lev=key>>(n-1)
            assert literal.admissible(lev,word,K);fields[kind,key]=(zl,zh,zl,zh);head_count+=1;continue
        word=literal.decode(key&((1<<m)-1),m);lev=key>>m;head=m//2;anchor=word;assert literal.admissible(lev,word,K)
        split='S1' if word[head]==1 else ('S2' if kind==1 else 'S3');labels=(split,'L' if kind==1 else 'R','M','R' if kind==1 else 'L')
        anchor_beta=cylinder(lev,anchor);ratio=None;displacement=[0,0]
        for label in labels:
            assert r-head>=0 and r-(len(word)-1-head)>=0
            assert word[:memory]==anchor[:memory] and word[-memory:]==anchor[-memory:] and sum(word)==sum(anchor)
            this=anchor_beta/cylinder(lev,word)
            if len(word)>m:
                if ratio is None:ratio=this
                assert ratio==this
            else:assert this==1
            dd=literal.DRIFT[label];displacement[0]+=dd[0];displacement[1]+=dd[1]
            if label=='L':head-=1
            elif label=='R':head+=1
            else:word,head=literal.literal_move(word,head,label)
            assert literal.admissible(lev,word,K)
        assert (word,head)==(anchor,m//2) and displacement==[0,0] and ratio>0
        ql,qh=enclose(point*ratio);fields[kind,key]=(zl,zh,ql,qh);cycle_count+=1
    assert len(descriptors)==len(coeff)==len(fields)
    del descriptors,coeff;gc.collect();event('cycle_identities',head_fields=head_count,edit_cycles=cycle_count)
    zero=g.mpfr(0);one=g.mpfr(1);Hlo=zero;Mlo=zero;Mhi=zero;min_density=None;windows_count=flows_count=0;max_relative_width=zero;extra=r-oldr
    flowmass={label:[zero,zero] for label in literal.LABELS};parity=[[zero,zero],[zero,zero]]
    for lev,w,c,blo,bhi in all_windows(n):
        oldword=w[extra:n-extra] if extra else w;oldlev=lev+sum(w[:extra]);oldcode=literal.code(oldword)
        vals=[];vl=zero;vh=zero
        for label in literal.legal_labels(lev,w,r,K):
            jl,jh=coarse[oldlev,oldcode,label]
            for kind,key,long in terms(lev,c,label,r):
                z=fields.get((kind,key))
                if z is not None:jl=down.add(jl,z[2 if long else 0]);jh=up.add(jh,z[3 if long else 1])
            assert jl>0 and jh>=jl,('nonpositive',lev,c,label,str(jl),str(jh))
            vals.append((label,jl,jh));vl=down.add(vl,jl);vh=up.add(vh,jh)
            min_density=jl if min_density is None or jl<min_density else min_density
            width=up.div(up.sub(jh,jl),jl)
            if width>max_relative_width:max_relative_width=width
        ml=down.mul(blo,vl);mh=up.mul(bhi,vh);Mlo=down.add(Mlo,ml);Mhi=up.add(Mhi,mh)
        pa=(lev+sum(w[:r]))&1;parity[pa][0]=down.add(parity[pa][0],ml);parity[pa][1]=up.add(parity[pa][1],mh)
        for label,jl,jh in vals:
            phi=min(one,up.div(jh,vl));entropy=down.mul(jl,-up.log(phi));Hlo=down.add(Hlo,down.mul(blo,entropy))
            fl=down.mul(blo,jl);fh=up.mul(bhi,jh);flowmass[label][0]=down.add(flowmass[label][0],fl);flowmass[label][1]=up.add(flowmass[label][1],fh)
        windows_count+=1;flows_count+=len(vals)
        if windows_count%20000==0:event('entropy',windows=windows_count,flows=flows_count)
    drift_intervals=[]
    for axis in (0,1):
        lo=zero;hi=zero
        for label,(fl,fh) in flowmass.items():
            if literal.DRIFT[label][axis]>0:lo=down.add(lo,fl);hi=up.add(hi,fh)
            else:lo=down.sub(lo,fh);hi=up.sub(hi,fl)
        assert lo<=0<=hi;drift_intervals.append([str(lo),str(hi)])
    assert down.sub(parity[0][0],parity[1][1])<=0<=up.sub(parity[0][1],parity[1][0])
    hlo=down.div(Hlo,Mhi);mu=down.exp(down.mul(g.mpfr(2),hlo));num,den=mu.as_integer_ratio();floor=Q((int(num)*10**12)//int(den),10**12)
    assert floor>=Q(61,5),('target not certified',str(floor))
    assert sha(a.fields)==fields_sha
    result=dict(status='PASS_EXPLICIT_CYCLES_EXACT_IDENTITIES_POSITIVITY_OUTWARD_ENTROPY',proof='Linear stationary-flow identity for an exact signed parent plus literal cycles; completed fine flow strictly positive',baseline_parent_status=manifest['status'],K=K,r=r,baseline_r=oldr,windows=windows_count,flows=flows_count,head_fields=head_count,edit_cycles=cycle_count,mu_lower_fraction=str(floor),mu_lower_decimal=f'{float(floor):.12f}',entropy_lower=str(hlo),mass_interval=[str(Mlo),str(Mhi)],minimum_density_lower=str(min_density),max_relative_density_interval_width=str(max_relative_width),drift_intervals=drift_intervals,precision=a.precision,mpfr_version=g.mpfr_version(),baseline_sha256=manifest['certificate_sha256'],baseline_manifest_sha256=sha(accepted),fields_sha256=fields_sha,verifier_sha256=sha(__file__),literal_sha256=sha(ROOT/'940_integer_literal_verifier/verify_sorted.py'),seconds=time.monotonic()-start)
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
