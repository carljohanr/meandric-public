"""Analytic Fisher projection onto reflected zero-mass head-cut directions."""
import numpy as np

class HeadCutProjection:
    def __init__(self,K,metric):
        self.K=K;self.sq=np.sqrt(metric);self.group=K.cut_orbit[K.component];self.weight=np.bincount(self.group,weights=metric[K.heads],minlength=K.cut_count);self.mass=np.bincount(self.group,weights=K.beta[K.heads],minlength=K.cut_count);self.den=float(np.sum(self.mass*self.mass/self.weight))
    def coefficients(self,x):
        total=np.bincount(self.group,weights=self.sq[self.K.heads]*x[self.K.heads],minlength=self.K.cut_count);lagrange=float(np.sum(self.mass*total/self.weight))/self.den;return (total-lagrange*self.mass)/self.weight
    def apply(self,x):
        z=self.coefficients(x);y=np.zeros_like(x);y[self.K.heads]=self.sq[self.K.heads]*z[self.group];return y
