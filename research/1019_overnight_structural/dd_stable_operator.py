"""Compensated forward and transpose basis products for ill-conditioned search coordinates."""
import numpy as np
from scipy.sparse import csr_matrix
from dd_column_kernel import DDColumnKernel
class DDStableKernel(DDColumnKernel):
 def __init__(self,*args,**kw):
  super().__init__(*args,**kw)
  self.column=self.B.tocsc();low=csr_matrix((self.low,self.B.indices,self.B.indptr),shape=self.B.shape).tocsc();assert np.array_equal(self.column.indices,low.indices) and np.array_equal(self.column.indptr,low.indptr)
  self.column_low=low.data.copy();self.column_pointer=self.column.indptr.astype(np.int32,copy=False);self.column_index=self.column.indices.astype(np.int32,copy=False)
 def apply(self,x):return self.apply_precise(x)
 def adjoint(self,g):
  g=(g+g[self.s.flowperm])/2;g=g-self.beta*(float(self.uniform@g)/self.mass_uniform);hg=g[self.heads];eg=g[self.edits]-self.E.T@self.td.adjoint(self.ta.adjoint(hg));eg=np.ascontiguousarray(np.bincount(self.edit_orbit,weights=eg,minlength=self.B.shape[0]));u=np.empty(self.B.shape[1]);self.multiply(len(u),self.column_pointer,self.column_index,self.column.data,self.column_low,eg,u);u-=self.size_gradient*float(self.size_gradient@u);z=np.bincount(self.cut_orbit[self.component],weights=hg,minlength=self.cut_count);return np.r_[u,z]
