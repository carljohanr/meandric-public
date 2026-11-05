"""Memory-light compensated selected-column evaluator, about106-bit arithmetic."""
from pathlib import Path
import ctypes,json,hashlib,sys
import numpy as np
from symmetric_head_kernel import SymmetricHeadKernel
class DDColumnKernel(SymmetricHeadKernel):
 def __init__(self,s,matrix_path,basis_path,echelon_path=None):
  super().__init__(s,matrix_path,basis_path);bp=Path(basis_path)
  with np.load(bp.with_name(bp.stem+'_metadata.npz')) as z:self.free=z['free'].copy()
  with np.load(bp.with_name(bp.stem+'_dd_low.npz')) as z:self.low=z['low'].copy();meta=json.loads(str(z['metadata']))
  assert self.low.shape==self.B.data.shape and meta['basis_sha256']==hashlib.file_digest(bp.open('rb'),'sha256').hexdigest()
  self.indptr=self.B.indptr.astype(np.int32,copy=False);self.indices=self.B.indices.astype(np.int32,copy=False)
  self.library=ctypes.CDLL(str(Path(__file__).parent/('double_double_csr.dylib' if sys.platform=='darwin' else 'double_double_csr.so')));self.multiply=self.library.dd_csr_mv
  I=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS');D=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');self.multiply.argtypes=[ctypes.c_int,I,I,D,D,D,D];self.multiply.restype=None
 def apply_precise(self,x):
  u=x[:self.B.shape[1]];u=np.ascontiguousarray(u-self.size_gradient*float(self.size_gradient@u));e=np.empty(self.B.shape[0]);self.multiply(len(e),self.indptr,self.indices,self.B.data,self.low,u,e);e=e[self.edit_orbit]
  h=self.reconstruct(e)+x[self.B.shape[1]:][self.cut_orbit[self.component]];v=np.zeros(self.shape[0]);v[self.edits]=e;v[self.heads]=h;v=(v+v[self.s.flowperm])/2;v-=self.uniform*(float(self.beta@v)/self.mass_uniform);return v
