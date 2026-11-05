"""Independent sparse literal-cycle replay over an accepted exact fine baseline.
Enumerates changed windows only; no full deeper model or stationarity matrix.
"""
import sys,json,pickle,hashlib,time,itertools,gc,argparse,resource
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'910_band_bound_push/vendor'));sys.path.insert(0,str(ROOT/'1019_overnight_structural'))
import gmpy2 as g
from signed_memory_cycle_replay import literal,terms,MemoryThreeLaw

def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
class Intervals:
 def __init__(self,precision):self.down=g.context(precision=precision,round=g.RoundDown);self.up=g.context(precision=precision,round=g.RoundUp);self.zero=g.mpfr(0)
 def enclose(self,q):
  if isinstance(q,Q):q=g.mpq(q.numerator,q.denominator)
  with self.down:lo=g.mpfr(q)
  with self.up:hi=g.mpfr(q)
  return lo,hi
 def add(self,a,b):return self.down.add(a[0],b[0]),self.up.add(a[1],b[1])
 def sub(self,a,b):return self.down.sub(a[0],b[1]),self.up.sub(a[1],b[0])
 def positive_times_signed(self,a,b):return self.down.mul(a[1] if b[0]<0 else a[0],b[0]),self.up.mul(a[0] if b[1]<0 else a[1],b[1])
 def xlogx(self,x):
  assert x[0]>0;lo=self.down.log(x[0]);hi=self.up.log(x[1]);return self.down.mul(x[1] if lo<0 else x[0],lo),self.up.mul(x[0] if hi<0 else x[1],hi)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--baseline-certificate',type=Path,required=True);ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--base-fields',type=Path,required=True);ap.add_argument('--increment',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--precision',type=int,default=160);ap.add_argument('--prior-increment',type=Path,action='append',default=[]);ap.add_argument('--inherited-certificate',type=Path);ap.add_argument('--parts',type=int,default=1);ap.add_argument('--part-range',type=int,nargs=2,default=None);ap.add_argument('--state-in',type=Path,default=None);ap.add_argument('--state-out',type=Path,default=None);a=ap.parse_args();assert not a.out.exists();NP=a.parts;PART=[0];start=time.monotonic();increment_digest=sha(a.increment);baseline_digest=sha(a.baseline_certificate);g.get_context().precision=a.precision;I=Intervals(a.precision);down,up=I.down,I.up;zero=(I.zero,I.zero)
 def event(stage,**kw):
  rec=dict(stage=stage,seconds=time.monotonic()-start,peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20);rec.update(kw);print(json.dumps(rec),flush=True)
 accepted=json.loads(a.baseline_certificate.read_text());rootaccepted=accepted;assert accepted['status']=='PASS_EXPLICIT_CYCLES_EXACT_IDENTITIES_POSITIVITY_OUTWARD_ENTROPY';assert sha(a.parent)==rootaccepted['baseline_sha256'] and sha(a.base_fields)==rootaccepted['fields_sha256'];assert g.mpfr(accepted['minimum_density_lower'])>0;rootaccepted=accepted;prior_hashes=[sha(p) for p in a.prior_increment];inherited_hash=None
 if a.prior_increment:
  assert a.inherited_certificate is not None;inherited_hash=sha(a.inherited_certificate);accepted=json.loads(a.inherited_certificate.read_text());assert accepted['status']=='PASS_SPARSE_LITERAL_CYCLES_INHERITED_POSITIVITY_OUTWARD_ENTROPY';assert accepted['baseline_certificate_sha256']==baseline_digest and accepted['increment_sha256']==prior_hashes[-1] and accepted.get('prior_increment_sha256',[])==prior_hashes[:-1]
 else:assert a.inherited_certificate is None
 r0=int(rootaccepted['r']);K=int(rootaccepted['K']);parent=pickle.load(a.parent.open('rb'));rp=int(parent['r']);params=parent['parameters'];assert len(params)==2 and parent['K']==K;lawup=parent['up'];kappa=parent['kappa'];law=MemoryThreeLaw(K,lawup,kappa,*params);PI={key:I.enclose(value) for key,value in law.forward.items()};initial={}
 def cylinder(lev,w):
  key=(lev,literal.code(w[:3]))
  if key not in initial:initial[key]=I.enclose(law.cylinder((lev,w[:3])))
  result=initial[key];vertices=[lev]
  for step in w[:3]:vertices.append(vertices[-1]+step)
  vertices=tuple(vertices)
  for step in w[3:]:
   nxt=vertices[-1]+step;result=I.positive_times_signed(result,PI[(*vertices,nxt)]);vertices=(*vertices[1:],nxt)
  return result
 coarse={};beta_cache={};D=I.enclose(parent['denominator']);np0=2*rp+1
 for key,N in zip(parent['keys'],parent['numerators']):
  key=int(key);wk=key>>3
  if wk not in beta_cache:beta_cache[wk]=cylinder(wk>>np0,literal.decode(wk&((1<<np0)-1),np0))
  blo,bhi=beta_cache[wk];nl,nh=I.enclose(N)
  if N>=0:coarse[key]=(down.div(down.div(nl,D[1]),bhi),up.div(up.div(nh,D[0]),blo))
  else:coarse[key]=(down.div(down.div(nl,D[0]),blo),up.div(up.div(nh,D[1]),bhi))
 del parent,beta_cache;gc.collect();event('accepted_coarse_density',flows=len(coarse))
 def read_fields(path,radius,need_support=False):
  n=2*radius+1;m=n-4;fields={};support=set()
  with np.load(path) as z:descriptors=z['descriptors'];coeff=z['coefficients']
  with np.load(path) as z:fm=json.loads(str(z['metadata']))
  if fm.get('family')=='template_trial_cycles':
   # Independent re-derivation of every instance: closed signed walks on literal configurations, sources fit the window.
   table={int(d['kind']):d for d in fm['templates']};direct={};unique=set();anchors=0
   for (kind,key),value in zip(descriptors,coeff):
    kind,key=int(kind),int(key)
    assert kind in table and np.isfinite(value) and (kind,key) not in unique;unique.add((kind,key));t=table[kind];tm=int(t['m']);th=int(t['h0']);lev=key>>tm;anchor=literal.decode(key&((1<<tm)-1),tm);assert literal.admissible(lev,anchor,K)
    point=Q.from_float(float(value));edges=[]
    if t['dtype']=='cycle':
     word=anchor;head=th
     for label in t['labels']:
      edges.append((word,head,label,1,word,head));word,head=literal.literal_move(word,head,label);edges[-1]=edges[-1][:4]+(word,head)
     assert (word,head)==(anchor,th)
    else:
     side,split=t['side'],t['split'];assert side in ('L','R') and split in ('S1','S2','S3')
     b=literal.literal_move(anchor,th,side);d=literal.literal_move(anchor,th,split);back='R' if side=='L' else 'L'
     c=(d[0],d[1]-1 if side=='L' else d[1]+1);assert literal.literal_move(c[0],c[1],back)==d and literal.literal_move(c[0],c[1],'M')==b
     edges=[(anchor,th,side,1,b[0],b[1]),(anchor,th,split,-1,d[0],d[1]),(c[0],c[1],'M',-1,b[0],b[1]),(c[0],c[1],back,1,d[0],d[1])]
    # circulation: net flow at every configuration is zero; protected moves; drift sums vanish
    net={};drift=[0,0]
    for (w0,h0_,lab,sg,w1,h1) in edges:
     assert literal.admissible(lev,w0,K) and literal.admissible(lev,w1,K) and w0[:3]==anchor[:3] and w0[-3:]==anchor[-3:] and w1[:3]==anchor[:3] and w1[-3:]==anchor[-3:] and sum(w0)==sum(w1)==sum(anchor)
     assert 0<h0_<len(w0)-1 and 0<h1<len(w1)-1
     net[(w0,h0_)]=net.get((w0,h0_),0)+sg;net[(w1,h1)]=net.get((w1,h1),0)-sg;dd=literal.DRIFT[lab];drift[0]+=sg*dd[0];drift[1]+=sg*dd[1]
    assert all(v==0 for v in net.values()) and drift==[0,0]
    weight=law.cylinder((lev,anchor))
    for (w0,h0_,lab,sg,w1,h1) in edges:
     left=radius-h0_;right=radius-(len(w0)-1-h0_);assert left>=0 and right>=0
     increment=I.enclose(point*sg*(weight/law.cylinder((lev,w0))))
     for exterior in itertools.product((1,-1),repeat=left+right):
      prefix=exterior[:left];full=prefix+w0+exterior[left:];startlev=lev-sum(prefix)
      if literal.admissible(startlev,full,K) and ((startlev<<n)|literal.code(full))%NP==PART[0]:
       packed=((((startlev<<n)|literal.code(full))<<3)|literal.LABELS.index(lab));direct[packed]=I.add(direct.get(packed,zero),increment);support.add(packed)
    anchors+=1
   return {'direct':direct},support
  if fm.get('family') in ('six_step_double_split_literal_cycles','eight_step_excursion_literal_cycles','mixed_literal_cycles'):
   patterns=(('L','S1','S3','R','M','M'),('L','S3','S3','R','M','M'),('R','S1','S2','L','M','M'),('R','S2','S2','L','M','M'))
   if fm.get('family')=='eight_step_excursion_literal_cycles':patterns=(('L', 'L', 'S1', 'R', 'R', 'R', 'M', 'L'), ('L', 'L', 'S3', 'R', 'R', 'R', 'M', 'L'), ('L', 'S1', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S1', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S1', 'S3', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S3', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'S3', 'R', 'R', 'M', 'M', 'L'), ('R', 'R', 'S1', 'L', 'L', 'L', 'M', 'R'), ('R', 'R', 'S2', 'L', 'L', 'L', 'M', 'R'), ('R', 'S1', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S1', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S1', 'S3', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S2', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'S3', 'L', 'L', 'M', 'M', 'R'))
   if fm.get('family')=='mixed_literal_cycles':patterns=patterns+(('L', 'L', 'S1', 'R', 'R', 'R', 'M', 'L'), ('L', 'L', 'S3', 'R', 'R', 'R', 'M', 'L'), ('L', 'S1', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S1', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S1', 'S3', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'R', 'R', 'R', 'M', 'L', 'L'), ('L', 'S3', 'S2', 'R', 'R', 'M', 'M', 'L'), ('L', 'S3', 'S3', 'R', 'R', 'M', 'M', 'L'), ('R', 'R', 'S1', 'L', 'L', 'L', 'M', 'R'), ('R', 'R', 'S2', 'L', 'L', 'L', 'M', 'R'), ('R', 'S1', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S1', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S1', 'S3', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'L', 'L', 'L', 'M', 'R', 'R'), ('R', 'S2', 'S2', 'L', 'L', 'M', 'M', 'R'), ('R', 'S2', 'S3', 'L', 'L', 'M', 'M', 'R'))
   anchorlen=n-6;direct={};unique=set()
   for (kind,key),value in zip(descriptors,coeff):
    kind,key=int(kind),int(key)
    if kind<3 and fm.get('family')=='mixed_literal_cycles':continue
    assert 3<=kind<3+len(patterns) and np.isfinite(value) and (kind,key) not in unique;unique.add((kind,key));lev=key>>anchorlen;anchor=literal.decode(key&((1<<anchorlen)-1),anchorlen);assert literal.admissible(lev,anchor,K);word=anchor;head=anchorlen//2;weight=law.cylinder((lev,anchor));displacement=[0,0];point=Q.from_float(float(value))
    for label in patterns[kind-3]:
     assert word[:3]==anchor[:3] and word[-3:]==anchor[-3:] and sum(word)==sum(anchor);ratio=weight/law.cylinder((lev,word));increment=I.enclose(point*ratio);left=radius-head;right=radius-(len(word)-1-head);assert left>=0 and right>=0
     for exterior in itertools.product((1,-1),repeat=left+right):
      prefix=exterior[:left];full=prefix+word+exterior[left:];startlev=lev-sum(prefix)
      if literal.admissible(startlev,full,K) and ((startlev<<n)|literal.code(full))%NP==PART[0]:
       packed=((((startlev<<n)|literal.code(full))<<3)|literal.LABELS.index(label));direct[packed]=I.add(direct.get(packed,zero),increment);support.add(packed)
     dd=literal.DRIFT[label];displacement[0]+=dd[0];displacement[1]+=dd[1];word,head=literal.literal_move(word,head,label);assert literal.admissible(lev,word,K)
    assert (word,head)==(anchor,anchorlen//2) and displacement==[0,0]
   if fm.get('family')!='mixed_literal_cycles':return {'direct':direct},support
   fields={'direct':direct}
  for (kind,key),value in zip(descriptors,coeff):
   kind,key=int(kind),int(key)
   if kind>=3 and fm.get('family')=='mixed_literal_cycles':continue
   assert kind in (0,1,2) and np.isfinite(value) and (kind,key) not in fields;point=Q.from_float(float(value));p=I.enclose(point)
   if kind==0:
    lev=key>>(n-1);word=literal.decode(key&((1<<(n-1))-1),n-1);assert literal.admissible(lev,word,K);fields[kind,key]=(*p,*p)
    if need_support:
     for step in (1,-1):
      full=word+(step,)
      if literal.admissible(lev,full,K) and ((lev<<n)|literal.code(full))%NP==PART[0]:support.add((((lev<<n)|literal.code(full))<<3)|0)
      full=(step,)+word;startlev=lev-step
      if literal.admissible(startlev,full,K) and ((startlev<<n)|literal.code(full))%NP==PART[0]:support.add((((startlev<<n)|literal.code(full))<<3)|1)
    continue
   lev=key>>m;anchor=literal.decode(key&((1<<m)-1),m);assert literal.admissible(lev,anchor,K);word=anchor;head=m//2;split='S1' if word[head]==1 else ('S2' if kind==1 else 'S3');labels=(split,'L' if kind==1 else 'R','M','R' if kind==1 else 'L');anchor_beta=law.cylinder((lev,anchor));ratio=None;displacement=[0,0]
   for label in labels:
    assert word[:3]==anchor[:3] and word[-3:]==anchor[-3:] and sum(word)==sum(anchor);this=anchor_beta/law.cylinder((lev,word))
    if len(word)>m:
     if ratio is None:ratio=this
     assert ratio==this
    else:assert this==1
    left=radius-head;right=radius-(len(word)-1-head);assert left>=0 and right>=0
    if need_support:
     for exterior in itertools.product((1,-1),repeat=left+right):
      prefix=exterior[:left];full=prefix+word+exterior[left:];startlev=lev-sum(prefix)
      if literal.admissible(startlev,full,K) and ((startlev<<n)|literal.code(full))%NP==PART[0]:support.add((((startlev<<n)|literal.code(full))<<3)|literal.LABELS.index(label))
    dd=literal.DRIFT[label];displacement[0]+=dd[0];displacement[1]+=dd[1]
    if label=='L':head-=1
    elif label=='R':head+=1
    else:word,head=literal.literal_move(word,head,label)
    assert literal.admissible(lev,word,K)
   assert (word,head)==(anchor,m//2) and displacement==[0,0] and ratio>0;q=I.enclose(point*ratio);fields[kind,key]=(*p,*q)
  return fields,support
 base_fields,_=read_fields(a.base_fields,r0);event('accepted_baseline_fields',fields=len(base_fields))
 with np.load(a.increment) as z:meta=json.loads(str(z['metadata']));savedkeys=z['changed_keys'].copy();saveddelta=z['density_increment'].copy()
 assert meta['status']=='NUMERICAL_ADAPTIVE_SPARSE_DEEP_CYCLE_NOT_CERTIFIED';r=int(meta['r']);n=2*r+1;assert r>r0 and meta['K']==K;source=Path(meta['source']);assert sha(source)==meta['source_sha256']
 cursor=source;last_r=r
 for prior in reversed(a.prior_increment):
  assert cursor.resolve()==prior.resolve()
  with np.load(prior) as z:pm=json.loads(str(z['metadata']))
  assert pm['status']=='NUMERICAL_ADAPTIVE_SPARSE_DEEP_CYCLE_NOT_CERTIFIED' and int(pm['r'])<=last_r and pm['K']==K;last_r=int(pm['r']);cursor=Path(pm['source']);assert sha(cursor)==pm['source_sha256']
 if a.prior_increment:
  with np.load(a.prior_increment[-1]) as z:assert int(accepted['r'])==int(json.loads(str(z['metadata']))['r'])
 with np.load(cursor) as z:sourcemeta=json.loads(str(z['metadata']))
 assert int(sourcemeta['r'])==r0 and 'memory_profile' not in sourcemeta and tuple(map(Q,sourcemeta['memory_parameters']))==params and {int(k):Q(v) for k,v in sourcemeta['up'].items()}==lawup and Q(sourcemeta['kappa'])==kappa
 total_support=0;total_windows=0;delta_H=zero;delta_M=zero;drift=[zero,zero];parity=zero;minimum=None;flows_checked=0;physical_proximity=0.;extra0=r0-rp;extraf=r-r0
 # EB1228 split-run support: the accumulators are exactly serialized between runs (mpfr as exact integer ratios), so a split run performs the same operations in the same order as a single run.
 P0,P1=(0,NP) if a.part_range is None else tuple(a.part_range)
 def dump_iv(x):return [list(map(str,u.as_integer_ratio())) for u in x]
 def load_iv(x):return tuple(g.mpfr(g.mpq(int(n),int(d))) for n,d in x)
 if a.state_in is not None:
  st=json.loads(a.state_in.read_text());assert st['increment_sha256']==increment_digest and st['parts']==NP and st['next_part']==P0 and st['precision']==a.precision and st['baseline_certificate_sha256']==baseline_digest
  total_support,total_windows,flows_checked,physical_proximity=st['total_support'],st['total_windows'],st['flows_checked'],st['physical_proximity'];delta_H,delta_M,parity=load_iv(st['delta_H']),load_iv(st['delta_M']),load_iv(st['parity']);drift=[load_iv(st['drift0']),load_iv(st['drift1'])];minimum=None if st['minimum'] is None else g.mpfr(g.mpq(int(st['minimum'][0]),int(st['minimum'][1])))
  for u,v in ((delta_H,st['delta_H']),(delta_M,st['delta_M']),(parity,st['parity'])):assert dump_iv(u)==v
 else:assert P0==0
 for PARTI in range(P0,P1):
  PART[0]=PARTI;event('part_start',part=PARTI)
  prior_fields=[];fields=None;support=None;gc.collect()  # EB1228: free the previous part before building the next (memory only)
  for path in a.prior_increment:
   with np.load(path) as z:pr=int(json.loads(str(z['metadata']))['r'])
   assert r0<pr<=r;pf,_=read_fields(path,pr);prior_fields.append((pr,pf))
  fields,support=read_fields(a.increment,r,True);pmask=((savedkeys>>3)%NP)==PARTI;assert np.array_equal(np.asarray(sorted(support),dtype=np.int64),savedkeys[pmask]);windows=sorted({key>>3 for key in support});saved={int(k):float(v) for k,v in zip(savedkeys[pmask],saveddelta[pmask])};total_support+=len(support);total_windows+=len(windows);event('independent_literal_support_and_identities',fields=len(fields),changed_flows=len(support),changed_windows=len(windows))
  def field_delta(lev,c,label,radius,authority):
   value=authority.get('direct',{}).get(((((lev<<(2*radius+1))|c)<<3)|literal.LABELS.index(label)),zero)
   for kind,key,long in terms(lev,c,label,radius):
    entry=authority.get((kind,key))
    if entry is not None:value=I.add(value,(entry[2 if long else 0],entry[3 if long else 1]))
   return value

  for wi,wk in enumerate(windows):
   lev=wk>>n;c=wk&((1<<n)-1);w=literal.decode(c,n);beta=cylinder(lev,w);oldlev=lev+sum(w[:extraf]);oldword=w[extraf:n-extraf];oldcode=literal.code(oldword);parentlev=oldlev+sum(oldword[:extra0]);parentword=oldword[extra0:len(oldword)-extra0] if extra0 else oldword;parentkey=((parentlev<<np0)|literal.code(parentword))<<3;v0=zero;v1=zero;dmwindow=zero
   for label in literal.legal_labels(lev,w,r,K):
    li=literal.LABELS.index(label);key=(wk<<3)|li;j0=I.add(coarse[parentkey|li],field_delta(oldlev,oldcode,label,r0,base_fields));
    for pr,pf in prior_fields:
     extra=r-pr;plev=lev+sum(w[:extra]);pword=w[extra:n-extra] if extra else w;j0=I.add(j0,field_delta(plev,literal.code(pword),label,pr,pf))
    delta=field_delta(lev,c,label,r,fields);j1=I.add(j0,delta);assert j0[0]>0 and j1[0]>0 and j1[1]>=j1[0],('nonpositive changed context',wk,label);minimum=j1[0] if minimum is None else min(minimum,j1[0]);v0=I.add(v0,j0);v1=I.add(v1,j1);flows_checked+=1
    if key not in support:assert delta==zero;continue
    f0=I.positive_times_signed(beta,j0);f1=I.positive_times_signed(beta,j1);df=I.positive_times_signed(beta,delta);delta_H=I.sub(delta_H,I.sub(I.xlogx(f1),I.xlogx(f0)));delta_M=I.add(delta_M,df);dmwindow=I.add(dmwindow,df)
    for axis in (0,1):drift[axis]=I.add(drift[axis],df) if literal.DRIFT[label][axis]>0 else I.sub(drift[axis],df)
    midpoint=float((delta[0]+delta[1])/2);physical_proximity+=float(beta[1])*abs(midpoint-saved[key])
   nu0=I.positive_times_signed(beta,v0);nu1=I.positive_times_signed(beta,v1);delta_H=I.add(delta_H,I.sub(I.xlogx(nu1),I.xlogx(nu0)));parity=I.add(parity,dmwindow) if (lev+sum(w[:r]))%2==0 else I.sub(parity,dmwindow)
   if wi%20000==0:event('sparse_outward_entropy',windows=wi,flows=flows_checked)
 if P1<NP:
  assert a.state_out is not None
  a.state_out.write_text(json.dumps(dict(increment_sha256=increment_digest,baseline_certificate_sha256=baseline_digest,parts=NP,next_part=P1,precision=a.precision,total_support=total_support,total_windows=total_windows,flows_checked=flows_checked,physical_proximity=physical_proximity,delta_H=dump_iv(delta_H),delta_M=dump_iv(delta_M),parity=dump_iv(parity),drift0=dump_iv(drift[0]),drift1=dump_iv(drift[1]),minimum=None if minimum is None else list(map(str,minimum.as_integer_ratio()))))+'\n');event('state_saved',next_part=P1);return
 assert physical_proximity<1e-8
 for interval in (*drift,parity):assert interval[0]<=0<=interval[1],('cycle invariant enclosure',interval)
 with down:h0=g.mpfr(accepted['entropy_lower']);M0lo=g.mpfr(accepted['mass_interval'][0])
 with up:M0hi=g.mpfr(accepted['mass_interval'][1])
 Hlo=down.add(down.mul(h0,M0lo),delta_H[0]);Mlo=down.add(M0lo,delta_M[0]);Mhi=up.add(M0hi,delta_M[1]);assert Hlo>0 and Mlo>0;hlo=down.div(Hlo,Mhi);mu=down.exp(down.mul(g.mpfr(2),hlo));num,den=mu.as_integer_ratio();floor=Q((int(num)*10**12)//int(den),10**12);assert abs(float(floor)-float(meta['mu']))<1e-7
 assert sha(a.parent)==rootaccepted['baseline_sha256'] and sha(a.base_fields)==rootaccepted['fields_sha256'] and sha(a.increment)==increment_digest and sha(a.baseline_certificate)==baseline_digest
 assert [sha(p) for p in a.prior_increment]==prior_hashes and (a.inherited_certificate is None or sha(a.inherited_certificate)==inherited_hash)
 result=dict(prior_increment_sha256=prior_hashes,inherited_certificate_sha256=inherited_hash,status='PASS_SPARSE_LITERAL_CYCLES_INHERITED_POSITIVITY_OUTWARD_ENTROPY',proof='Accepted exact fine baseline plus independently verified literal closed-cycle fields; untouched positivity inherited, every affected context positive, entropy/mass change outward.',K=K,r=r,baseline_r=r0,changed_flows=total_support,changed_windows=total_windows,flows_in_changed_windows_checked=flows_checked,full_deep_replay=False,mu_lower_decimal=f'{float(floor):.12f}',mu_lower_fraction=str(floor),entropy_lower=str(hlo),mass_interval=[str(Mlo),str(Mhi)],minimum_checked_density_lower=str(minimum),increment_drift_intervals=[[str(x),str(y)] for x,y in drift],increment_parity_interval=list(map(str,parity)),physical_l1_increment_rounding_difference=physical_proximity,precision=a.precision,baseline_certificate_sha256=sha(a.baseline_certificate),parent_sha256=sha(a.parent),base_fields_sha256=sha(a.base_fields),increment_sha256=sha(a.increment),verifier_sha256=sha(__file__),seconds=time.monotonic()-start);a.out.write_text(json.dumps(result,indent=2)+'\n');event('finished',**result)
if __name__=='__main__':main()
