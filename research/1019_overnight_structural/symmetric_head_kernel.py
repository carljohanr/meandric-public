"""Reflection-first complete stationary kernel with exact edit-basis provenance."""
import pickle
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz
from scipy.sparse.linalg import LinearOperator
from head_compatibility import compatibility

class SymmetricHeadKernel:
    def __init__(self,s,matrix_path,basis_path):
        recipe=pickle.load(Path(matrix_path).open('rb'));recipe.pop('integer_rows',None);assert recipe['metadata']['reflect'] and not recipe['metadata']['include_size'];assert np.array_equal(recipe['keys'],s.keys)
        self.s=s;self.B=load_npz(basis_path);C,self.heads,self.edits,self.reconstruct,self.td,self.ta=compatibility(s);assert np.array_equal(recipe['full_edits'],self.edits);self.edit_orbit=recipe['edit_orbit'];assert self.B.shape[0]==len(recipe['edits']);self.E=s.M[:,self.edits].tocsr()
        self.component=np.full(len(self.heads),-1,dtype=int);self.component[self.ta.roots]=np.arange(len(self.ta.roots))
        for v in self.ta.levels:self.component[v]=self.component[self.ta.parent[v]]
        head_index=np.full(len(s.keys),-1,dtype=int);head_index[self.heads]=np.arange(len(self.heads));partner=self.component[head_index[s.flowperm[self.heads[self.ta.roots]]]];representatives=np.minimum(np.arange(len(partner)),partner);unique,self.cut_orbit=np.unique(representatives,return_inverse=True);self.cut_count=len(unique)
        self.beta=s.bw[s.vw];self.uniform=((s.keys&7)<2).astype(float);self.mass_uniform=float(self.beta@self.uniform);edit_size=np.bincount(self.edit_orbit,weights=self.beta[self.edits]*s.db[self.edits]);self.size_gradient=self.B.T@edit_size;self.size_gradient/=max(np.linalg.norm(self.size_gradient),1e-300)
        exact_size_path=Path(basis_path).with_name(Path(basis_path).stem+'_exact_size.npz')
        if exact_size_path.exists():
            import hashlib,json
            with np.load(exact_size_path) as exact_size:
                authority=json.loads(str(exact_size['metadata']));assert authority['basis_sha256']==hashlib.file_digest(Path(basis_path).open('rb'),'sha256').hexdigest();self.size_gradient=exact_size['gradient'].copy();self.size_gradient_authority=authority
        self.shape=(len(s.keys),self.B.shape[1]+self.cut_count);self.operator=LinearOperator(self.shape,matvec=self.apply,rmatvec=self.adjoint,dtype=float)
    def apply(self,x):
        u=x[:self.B.shape[1]];u=u-self.size_gradient*float(self.size_gradient@u);e=(self.B@u)[self.edit_orbit];h=self.reconstruct(e)+x[self.B.shape[1]:][self.cut_orbit[self.component]];v=np.zeros(self.shape[0]);v[self.edits]=e;v[self.heads]=h;v=(v+v[self.s.flowperm])/2;v-=self.uniform*(float(self.beta@v)/self.mass_uniform);return v
    def adjoint(self,g):
        g=(g+g[self.s.flowperm])/2;g=g-self.beta*(float(self.uniform@g)/self.mass_uniform);hg=g[self.heads];eg=g[self.edits]-self.E.T@self.td.adjoint(self.ta.adjoint(hg));u=self.B.T@np.bincount(self.edit_orbit,weights=eg,minlength=self.B.shape[0]);u-=self.size_gradient*float(self.size_gradient@u);z=np.bincount(self.cut_orbit[self.component],weights=hg,minlength=self.cut_count);return np.r_[u,z]
