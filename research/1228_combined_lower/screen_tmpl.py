"""Stream template candidates against the record baseline; keep best-scored directions per group (numerical only)."""
from tmpl import *
import argparse,hashlib,resource,time
ap=argparse.ArgumentParser()
ap.add_argument('--source',type=Path,required=True);ap.add_argument('--r',type=int,required=True)
ap.add_argument('--groups',default='c4,e8,dia');ap.add_argument('--keep',default='');ap.add_argument('--draws',type=int,default=300000,help='random draws per template')
ap.add_argument('--batch',type=int,default=512);ap.add_argument('--law-fraction',type=float,default=.7);ap.add_argument('--seed',type=int,default=1211)
ap.add_argument('--children-of',type=Path,default=None,help='record npz: add one-letter refinements of its production four-step descriptors (c4 group)')
ap.add_argument('--out',type=Path,required=True)
a=ap.parse_args();assert not a.out.exists();start=time.monotonic()
groups=tuple(a.groups.split(','));keep={g:16384 for g in groups}
for item in filter(None,a.keep.split(',')):g,v=item.split('=');keep[g]=int(v)
base=SparseFineBaseline(a.source,a.r);G=TemplateGenerator(base.K,a.r,base.law,base.parameters,groups=groups);rng=np.random.default_rng(a.seed)
def group_of(t):return t.name.split(':')[0]
heaps={g:[] for g in groups};seen={g:0 for g in groups};legal={g:0 for g in groups};stats={t.name:[0,0] for t in G.templates}
def emit(stage,**kw):
    rec=dict(stage=stage,seconds=round(time.monotonic()-start,1),peak_mib=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,1));rec.update(kw);print(json.dumps(rec),flush=True)
emit('baseline',mu=base.initial_mu,templates=len(G.templates))
def candidates():
    done=set()
    if a.children_of is not None:
        with np.load(a.children_of) as z:desc=z['descriptors']
        byname={t.name:t for t in G.templates};m0=2*a.r-3;h=m0//2;K=base.K
        for kind,key in desc:
            kind,key=int(kind),int(key)
            if kind not in (1,2):continue
            lev=key>>m0;word=literal.decode(key&((1<<m0)-1),m0)
            if not literal.admissible(lev,word,K):continue
            for ell in (1,-1):
                if kind==1:
                    nw=(ell,)+word;nl=lev-ell;name='c4:S1 L M R' if word[h]==1 else 'c4:S2 L M R'
                else:
                    nw=word+(ell,);nl=lev;name='c4:S1 R M L' if word[h]==1 else 'c4:S3 R M L'
                if name not in byname or not literal.admissible(nl,nw,K):continue
                t=byname[name];d=(t.kind,(nl<<t.m)|literal.code(nw))
                if d in done:continue
                done.add(d);e=G.entries(d)
                if e:yield d,e
        emit('children_generated',count=len(done))
    for t in G.templates:
        target=a.draws
        for offset in range(0,target,4096):
            keys=sample_window_keys(base,t.m,min(4096,target-offset),rng,a.law_fraction)
            for key in keys:
                d=(t.kind,int(key))
                if d in done:continue
                done.add(d);stats[t.name][0]+=1;e=G.entries(d)
                if e:stats[t.name][1]+=1;yield d,e
it=iter(candidates());batches=0;names={t.kind:t.name for t in G.templates}
while True:
    cand=list(itertools.islice(it,a.batch))
    if not cand:break
    keys=array('q');cols=array('i');vals=array('d');C=len(cand)
    for col,(_,e) in enumerate(cand):
        for k,v in e:keys.append(k);cols.append(col);vals.append(v)
    keys=np.asarray(keys);cols=np.asarray(cols);v=np.asarray(vals);beta,j,J=base.values(keys)
    g=np.bincount(cols,weights=beta*v*(-np.log(j/J)-base.h),minlength=C);curv=np.bincount(cols,weights=beta*v*v/j,minlength=C)
    pairs=(keys>>3)*C+cols;unique,first,group=np.unique(pairs,return_index=True,return_inverse=True);sums=np.bincount(group,weights=v)
    curv-=np.bincount(unique%C,weights=beta[first]*sums*sums/J[first],minlength=C);curv=np.maximum(curv,1e-300)
    with np.errstate(divide='ignore'):
        ratio=np.where(v>0,j/np.abs(v),np.inf);fneg=np.full(C,np.inf);np.minimum.at(fneg,cols,ratio)
        ratio=np.where(v<0,j/np.abs(v),np.inf);fpos=np.full(C,np.inf);np.minimum.at(fpos,cols,ratio)
    step=g/curv;step=np.clip(step,-.8*fneg,.8*fpos);score=np.maximum(0,g*step-.5*curv*step*step)
    for col,(d,_) in enumerate(cand):
        gname=names[d[0]].split(':')[0];seen[gname]+=1;item=(float(score[col]),int(d[0]),int(d[1]),float(g[col]),float(curv[col]),float(min(fneg[col],fpos[col])));hp=heaps[gname]
        if len(hp)<keep[gname]:heapq.heappush(hp,item)
        elif item[0]>hp[0][0]:heapq.heapreplace(hp,item)
    batches+=1
    if batches%200==0:emit('progress',seen=seen,retained={k:len(x) for k,x in heaps.items()})
chosen=[];per={}
for gname,hp in heaps.items():
    s=sorted(hp,reverse=True);chosen+=s;per[gname]=dict(retained=len(s),top_score=s[0][0] if s else 0,median_score=s[len(s)//2][0] if s else 0,sum_score=sum(x[0] for x in s))
chosen.sort(reverse=True)
res=dict(family='template_trial_cycles',status='NUMERICAL_STREAMED_ADAPTIVE_DEEP_CYCLE_SELECTION',r=a.r,K=base.K,baseline_mu=base.initial_mu,source=str(a.source),source_sha256=hashlib.file_digest(a.source.open('rb'),'sha256').hexdigest(),groups=list(groups),keep=keep,draws_per_template=a.draws,law_fraction=a.law_fraction,children_of=str(a.children_of) if a.children_of else None,per_group=per,template_stats=stats,templates=[t.as_dict() for t in G.templates],descriptors=[[x[1],x[2]] for x in chosen],scores=[x[0] for x in chosen],seconds=time.monotonic()-start,scope='Numerical selection only.')
a.out.write_text(json.dumps(res)+'\n');emit('finished',per_group=per,seen=seen,template_stats=stats)
