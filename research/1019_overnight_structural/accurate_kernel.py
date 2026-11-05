"""High-precision final edit reconstruction from the exact triangular recipe.

The iterative operator stays sparse binary64. Final directions are reconstructed
at 50 decimal digits before conversion, to avoid amplified basis rounding.
"""
import pickle,json,time
from pathlib import Path
from decimal import Decimal,localcontext
import numpy as np
from symmetric_head_kernel import SymmetricHeadKernel

class AccurateSymmetricHeadKernel(SymmetricHeadKernel):
    def __init__(self,s,matrix_path,basis_path,echelon_path):
        super().__init__(s,matrix_path,basis_path)
        recipe=pickle.load(Path(echelon_path).open('rb'));self.pivots=recipe['pivots'];priority=recipe['priority'];self.back_order=sorted(self.pivots,key=lambda c:priority[c],reverse=True)
        with np.load(Path(basis_path).with_name(Path(basis_path).stem+'_metadata.npz')) as z:
            self.free=z['free'].copy()
            with localcontext() as ctx:
                ctx.prec=50;self.exact_scales=[Decimal(str(n))/Decimal(str(d)) for n,d in zip(z['scale_numerators'],z['scale_denominators'])]
    def apply_precise(self,x):
        u=x[:self.B.shape[1]];u=u-self.size_gradient*float(self.size_gradient@u)
        with localcontext() as ctx:
            ctx.prec=50;e=[Decimal(0)]*self.B.shape[0]
            for c,value,scale in zip(self.free,u,self.exact_scales):e[c]=Decimal.from_float(float(value))/scale
            for c in self.back_order:
                row=self.pivots[c];e[c]=-sum((Decimal(v)*e[j] for j,v in row.items() if j!=c),Decimal(0))/Decimal(row[c])
            e=np.asarray([float(v) for v in e])[self.edit_orbit]
        h=self.reconstruct(e)+x[self.B.shape[1]:][self.cut_orbit[self.component]];v=np.zeros(self.shape[0]);v[self.edits]=e;v[self.heads]=h;v=(v+v[self.s.flowperm])/2;v-=self.uniform*(float(self.beta@v)/self.mass_uniform);return v
