"""Matched radius-seven entropy screen of the complete sparse r5 kernel."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'1017_transport_flows'))
from search import *
from fractions import Fraction
from dd_column_kernel import DDColumnKernel as SymmetricHeadKernel
from head_cut_projection import HeadCutProjection
from scipy.sparse.linalg import LinearOperator,lsmr,cg,splu
from scipy.sparse import coo_matrix,diags
from scipy.linalg import cholesky,solve_triangular
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--warm-increment',type=Path);ap.add_argument('--warm-basis-prefix',type=Path);ap.add_argument('--positive-parent',type=Path,required=True);ap.add_argument('--matrix-prefix',type=Path,required=True);ap.add_argument('--basis-prefix',type=Path,required=True);ap.add_argument('--coarse-r',type=int,default=6);ap.add_argument('--block',type=int,default=1024);ap.add_argument('--cg-iterations',type=int,default=500);ap.add_argument('--rounds',type=int,default=4);ap.add_argument('--head-schur',action='store_true');ap.add_argument('--coupled-preconditioner',action='store_true');ap.add_argument('--accurate-operator',action='store_true');ap.add_argument('--remove-parameter-gauges',action='store_true');ap.add_argument('--trust-floor',type=float,default=0.0);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.with_suffix('.json').exists()
if a.accurate_operator:
 from dd_stable_operator import DDStableKernel as SymmetricHeadKernel
p=Path(__file__).parent
with np.load(a.witness,allow_pickle=False) as source:
    metadata=json.loads(str(source['metadata']));up={int(i):Fraction(v) for i,v in metadata['up'].items()};kappa=Fraction(metadata['kappa']);order=int(metadata.get('memory_order',2));eta=Fraction(metadata['memory_eta']);Kheight=int(metadata['K']);fine_r=int(metadata['r']);physical=source['physical'].copy();source_keys=source['keys'].copy()
from local_ratio_builder import build_local
parameters=tuple(Fraction(v) for v in metadata.get('memory_parameters',[metadata.get('base_eta','1/8'),str(eta)] if order==3 else [str(eta)]))
def make(r):return build_local(Kheight,r,up,kappa,parameters)
from local_cycle_family import geometry_witness
s=geometry_witness(a.witness,fine_r) if fine_r>a.coarse_r else make(fine_r);assert np.array_equal(source_keys,s.keys);beta=s.bw[s.vw];j=physical/beta;assert np.all(j>0)
cs=make(a.coarse_r);K=SymmetricHeadKernel(cs,a.matrix_prefix.with_name(a.matrix_prefix.name+'_matrix.pkl'),a.basis_prefix.with_suffix('.npz'),a.matrix_prefix.with_name(a.matrix_prefix.name+'_echelon.pkl'));mapping=restriction(s.keys,s.r,cs.keys,a.coarse_r);initial_j=j.copy();parent_geometry=geometry_witness(a.positive_parent,a.coarse_r);assert np.array_equal(parent_geometry.keys,cs.keys);parent_base=parent_geometry.base.copy();del parent_geometry;theta=np.zeros(K.shape[1]);records=[];rng=np.random.default_rng(20261012)

def score():
    f=beta*j;nu=np.bincount(s.vw,weights=f,minlength=len(s.bw));logs=np.log(f/nu[s.vw]);mass=f.sum();h=-f@logs/mass;return float(np.exp(2*h)),h,f,nu,logs,mass
if a.warm_increment:
    assert a.warm_basis_prefix
    with np.load(a.warm_increment) as warm:old_effective=warm['effective_edit_coefficients'].copy();old_theta=warm['coefficients'].copy()
    with np.load(a.warm_basis_prefix.with_name(a.warm_basis_prefix.name+'_metadata.npz')) as warm_basis:old_free=warm_basis['free'].copy()
    position={int(c):i for i,c in enumerate(K.free)}
    for c,v in zip(old_free,old_effective):theta[position[int(c)]]=v
    theta[K.B.shape[1]:]=old_theta[len(old_free):];j=initial_j+K.apply_precise(theta)[mapping];assert np.all(j>0)
initial=score()[0];event('sparse_kernel_screen_start',initial_mu=initial,coarse_flows=K.shape[0],parameters=K.shape[1])
with np.load(a.basis_prefix.with_name(a.basis_prefix.name+'_metadata.npz')) as basis_metadata:representative_edits=basis_metadata['edits'].copy()
for turn in range(a.rounds):
    old,h,f,nu,logs,mass=score();g=np.bincount(mapping,weights=beta*(-logs-h),minlength=len(cs.keys));metric=np.bincount(mapping,weights=beta/j+a.trust_floor/(len(j)*j*j),minlength=len(cs.keys));sq=np.sqrt(metric);head_projection=HeadCutProjection(K,metric);edit_count=K.B.shape[1];target=g/sq;head_fixed=head_projection.apply(target);head_fixed_coefficients=head_projection.coefficients(target);residual_target=target-head_fixed;diagonal=np.zeros(edit_count)
    def edit_forward(u):
        v=sq*K.apply(np.r_[u,np.zeros(K.cut_count)]);return v-head_projection.apply(v)
    def edit_adjoint(v):return K.adjoint(sq*(v-head_projection.apply(v)))[:edit_count]
    before=time.monotonic();edit_weight=np.bincount(K.edit_orbit,weights=metric[K.edits],minlength=K.B.shape[0]);norms=np.sqrt(K.B.multiply(K.B).T@edit_weight);scale=1/np.maximum(norms,1e-100);column_order=np.argsort(cs.keys[representative_edits[K.free]]);blocks=[]
    for offset in range(0,edit_count,a.block):
        ids=column_order[offset:offset+a.block];bb=K.B[:,ids];gram=(bb.T@bb.multiply(edit_weight[:,None])).toarray();gram*=scale[ids,None]*scale[None,ids]
        for ridge in (1e-14,1e-12,1e-10,1e-8,1e-6):
            try:factor=cholesky(gram+ridge*np.eye(len(ids)),lower=True,check_finite=False);break
            except np.linalg.LinAlgError:continue
        blocks.append((ids,factor))
    def transform(x):
        out=np.empty(edit_count)
        for ids,factor in blocks:out[ids]=scale[ids]*solve_triangular(factor,x[ids],lower=True,trans='T',check_finite=False)
        return out
    def transform_adjoint(x):
        out=np.empty(edit_count)
        for ids,factor in blocks:out[ids]=solve_triangular(factor,scale[ids]*x[ids],lower=True,check_finite=False)
        return out
    edit_transform=transform;edit_transform_adjoint=transform_adjoint;head_scale=1/np.sqrt(head_projection.weight)
    def transform_full(x):
        u=edit_transform(x[:edit_count]);z=head_scale*x[edit_count:]
        if a.head_schur:z-=head_projection.coefficients(sq*K.apply(np.r_[u,np.zeros(K.cut_count)]))
        return np.r_[u,z]
    def adjoint_full(x):
        ge=x[:edit_count].copy();gh=x[edit_count:]
        if a.head_schur:
            mw=head_projection.mass/head_projection.weight;t=gh/head_projection.weight-mw*float(mw@gh)/head_projection.den;v=np.zeros(K.shape[0]);v[K.heads]=metric[K.heads]*t[head_projection.group];ge-=K.adjoint(v)[:edit_count]
        return np.r_[edit_transform_adjoint(ge),head_scale*gh]
    xx=rng.normal(size=K.shape[1]);yy=rng.normal(size=K.shape[1]);aa=float(transform_full(xx)@yy);bb=float(xx@adjoint_full(yy));adjerr=abs(aa-bb)/max(1,abs(aa),abs(bb));assert adjerr<1e-9;event('block_transform_control',relative_adjoint_error=adjerr,blocks=len(blocks),seconds_prepare=time.monotonic()-before)
    first=np.searchsorted(s.vw,np.arange(len(s.bw)));window_mapping=cs.vw[mapping[first]];vertex_weight=np.bincount(window_mapping,weights=s.bw*s.bw/nu,minlength=len(cs.bw))
    def hessian(v):return metric*v-vertex_weight[cs.vw]*np.bincount(cs.vw,weights=v,minlength=len(cs.bw))[cs.vw]
    group=K.cut_orbit[K.component];wh=cs.vw[K.heads];incidence=coo_matrix((np.ones(len(K.heads)),(wh,group)),shape=(len(cs.bw),K.cut_count)).tocsr();head_matrix=diags(head_projection.weight)-incidence.T@incidence.multiply(vertex_weight[:,None]);hs=1/np.sqrt(head_projection.weight);scaled_head=head_matrix.multiply(hs[:,None]).multiply(hs[None,:]).tocsc();head_factor=splu(scaled_head+1e-10*diags(np.ones(K.cut_count)),permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0);vm=hs*head_factor.solve(hs*head_projection.mass);den=float(head_projection.mass@vm)
    def head_inverse(x):
        z=hs*head_factor.solve(hs*x);return z-vm*float(head_projection.mass@z)/den
    def precondition(x):
        zh=head_inverse(x[edit_count:])
        ge=x[:edit_count]
        if a.coupled_preconditioner:
            ge=ge-K.adjoint(hessian(K.apply(np.r_[np.zeros(edit_count),zh])))[:edit_count]
        ze=edit_transform(edit_transform_adjoint(ge))
        if a.coupled_preconditioner:
            zh=zh-head_inverse(K.adjoint(hessian(K.apply(np.r_[ze,np.zeros(K.cut_count)])))[edit_count:])
        return np.r_[ze,zh]
    xx=rng.normal(size=K.shape[1]);yy=rng.normal(size=K.shape[1]);px=precondition(xx);py=precondition(yy);symerr=abs(float(xx@py-yy@px))/max(1,abs(float(xx@py)),abs(float(yy@px)));assert symerr<1e-7 and xx@px>0;event('conditional_head_preconditioner',symmetry_error=symerr,head_factor_nnz=int(head_factor.L.nnz+head_factor.U.nnz),seconds_prepare=time.monotonic()-before)
    def parameter_project(x):
        out=x.copy();q=K.size_gradient;qn=float(q@q)
        if qn:out[:edit_count]-=q*float(q@out[:edit_count])/qn
        out[edit_count:]-=float(np.mean(out[edit_count:]))
        return out
    def raw_hessian(x):return K.adjoint(hessian(K.apply(x)))
    if a.remove_parameter_gauges:
        def projected_hessian(x):
            px=parameter_project(x);return parameter_project(raw_hessian(px))+(x-px)
        def projected_precondition(x):
            px=parameter_project(x);return parameter_project(precondition(px))+(x-px)
        H=LinearOperator((K.shape[1],K.shape[1]),matvec=projected_hessian,dtype=float);M=LinearOperator(H.shape,matvec=projected_precondition,dtype=float);rhs=parameter_project(K.adjoint(g))
        qtest=np.r_[K.size_gradient,np.zeros(K.cut_count)];htest=np.r_[np.zeros(edit_count),np.ones(K.cut_count)];event('parameter_gauge_control',edit_null_max=float(np.max(abs(K.apply(qtest)))),head_null_max=float(np.max(abs(K.apply(htest)))),projector_error=float(np.max(abs(parameter_project(parameter_project(xx))-parameter_project(xx)))))
    else:H=LinearOperator((K.shape[1],K.shape[1]),matvec=raw_hessian,dtype=float);M=LinearOperator(H.shape,matvec=precondition,dtype=float);rhs=K.adjoint(g)
    counter=[0]
    def callback(x):counter[0]+=1
    solution,info=cg(H,rhs,M=M,rtol=1e-9,atol=0,maxiter=a.cg_iterations,callback=callback);dt=solution;keep=np.ones(edit_count,dtype=bool);sol=(None,int(info),counter[0]);u=dt[:edit_count];event('conditional_newton',cg_info=int(info),iterations=counter[0],relative_linear_residual=float(np.linalg.norm(H@solution-rhs)/np.linalg.norm(rhs)),seconds_solve=time.monotonic()-before);coarse=K.apply_precise(dt);fine=coarse[mapping];negative=fine<0;maximum=float(.98*np.min(-j[negative]/fine[negative])) if negative.any() else 1.;df=beta*fine;dn=np.bincount(s.vw,weights=df,minlength=len(nu));linear=float(df@(-logs-h))
    def gain(t):
        return (t*linear-f@phi(t*fine/j)+nu@phi(t*dn/nu))/(mass+t*df.sum())
    bracket=min(1.,maximum)
    for _ in range(20):
        if bracket>=maximum or gain(bracket)<=gain(bracket/2):break
        bracket=min(maximum,2*bracket)
    opt=minimize_scalar(lambda t:-gain(t),bounds=(0,bracket),method='bounded',options=dict(xatol=1e-10,maxiter=70));step=float(opt.x);inc=gain(step);raw=float(max(np.max(abs(cs.M@(step*coarse))),np.max(abs(cs.M@K.apply_precise(theta+step*dt)))));bc=cs.bw[cs.vw];other=float(max(abs((bc*cs.da)@(step*coarse)),abs((bc*cs.db)@(step*coarse)),abs(bc@(step*coarse))));trial=initial_j+K.apply_precise(theta+step*dt)[mapping]
    if s.r==cs.r:raw=max(raw,float(np.max(abs(s.M@trial))))
    accepted=bool(inc>1e-13 and raw<1e-8 and other<1e-10 and np.all(trial>0))
    for back in range(12):
        if accepted or inc<=1e-13:break
        step*=.5;inc=gain(step);raw=float(max(np.max(abs(cs.M@(step*coarse))),np.max(abs(cs.M@K.apply_precise(theta+step*dt)))));other=float(max(abs((bc*cs.da)@(step*coarse)),abs((bc*cs.db)@(step*coarse)),abs(bc@(step*coarse))));trial=initial_j+K.apply_precise(theta+step*dt)[mapping]
        if s.r==cs.r:raw=max(raw,float(np.max(abs(s.M@trial))))
        accepted=bool(inc>1e-13 and raw<1e-8 and other<1e-10 and np.all(trial>0))
    if accepted:
        theta+=step*dt;checkpoint_direction=K.apply_precise(theta);j=initial_j+checkpoint_direction[mapping];checkpoint_meta=dict(metadata);checkpoint_meta.update(r=a.coarse_r,status='NUMERICAL_SIGNED_PARENT_INCREMENT_CHECKPOINT');effective_checkpoint=theta[:K.B.shape[1]]-K.size_gradient*float(K.size_gradient@theta[:K.B.shape[1]]);np.savez(a.out.with_name(a.out.name+'_parent_checkpoint.npz'),keys=cs.keys,physical=cs.bw[cs.vw]*(parent_base+checkpoint_direction),effective_edit_coefficients=effective_checkpoint,coefficients=theta,metadata=np.array(json.dumps(checkpoint_meta)))
    rec=dict(round=turn,old_mu=old,mu=score()[0],gain_h=inc,step=step,maximum=maximum,bracket=bracket,accepted=accepted,density_residual=raw,drift_mass_residual=other,lsmr_iterations=int(sol[2]),lsmr_stop=int(sol[1]),solve_seconds=time.monotonic()-before,parameters_retained=int(keep.sum()));records.append(rec);event('sparse_kernel_proposal',**rec)
    if not accepted or inc<1e-10:break
aggregate=K.apply_precise(theta);result=dict(K=Kheight,r=fine_r,kappa=str(kappa),up={str(i):str(v) for i,v in up.items()},memory_parameters=list(map(str,parameters)),memory_order=order,base_eta=metadata.get('base_eta'),memory_eta=str(eta),source_witness=str(a.witness),coarse_r=a.coarse_r,block=a.block,head_schur=a.head_schur,coupled_preconditioner=a.coupled_preconditioner,accurate_operator=a.accurate_operator,trust_floor=a.trust_floor,preconditioner='conditional_head_sparse_inverse_and_edit_blocks',status='COMPLETE_SPARSE_KERNEL_NUMERICAL_ENTROPY_SCREEN',initial_mu=initial,final_mu=score()[0],records=records,seconds=time.monotonic()-START,aggregate_coarse_density_residual=float(np.max(abs(cs.M@aggregate))),actual_density_residual=float(np.max(abs(s.M@j))) if s.r==cs.r else None,scope='Sparse exact-space provenance with floating coefficients, original coarse-row checks, no full fine-radius replay or certificate.')
effective=theta[:K.B.shape[1]]-K.size_gradient*float(K.size_gradient@theta[:K.B.shape[1]]);parent=parent_base+aggregate;parent_meta=dict(metadata);parent_meta.update(r=a.coarse_r,status='NUMERICAL_SIGNED_COARSE_PARENT_EXACT_BASE_PLUS_SELECTED_INCREMENT');np.savez(a.out.with_name(a.out.name+'_parent.npz'),keys=cs.keys,physical=cs.bw[cs.vw]*parent,metadata=np.array(json.dumps(parent_meta)));np.savez(a.out.with_suffix('.npz'),effective_edit_coefficients=effective,keys=s.keys,physical=score()[2],coefficients=theta,coarse_keys=cs.keys,coarse_direction=aggregate,metadata=np.array(json.dumps(result)));a.out.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');event('sparse_kernel_screen_finished',mu=result['final_mu'])
