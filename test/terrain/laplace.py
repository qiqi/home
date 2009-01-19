import numpy
from numpy import zeros, ones, logical_or, log2

def laplace(phi):
   # laplace operator
   result = zeros(phi.shape, dtype=float)
   result[1:,:] += phi[:-1,:]
   result[:-1,:] += phi[1:,:]
   result[:,1:] += phi[:,:-1]
   result[:,:-1] += phi[:,1:]
   result[:,:] -= 4.0 * phi[:,:]
   return result

def resid(phi, source, bcmask):
   res = source - laplace(phi)
   res[bcmask] = 0.0
   return res

def collcect(res, mask):
   assert res.shape == mask.shape
   nx, ny = res.shape
   assert nx % 2 == ny % 2 == 0
   res_c = res[0::2,0::2] + res[1::2,0::2] + res[0::2,1::2] + res[1::2,1::2]
   mask_c = logical_or(logical_or(mask[0::2,0::2], mask[1::2,0::2]), \
                       logical_or(mask[0::2,1::2], mask[1::2,1::2]))
   return res_c, mask_c

def distribute(phi_c):
   nx, ny = phi_c.shape
   phi = zeros([2*nx, 2*ny])
   phi[0::2,0::2] = phi_c
   phi[1::2,0::2] = phi_c
   phi[0::2,1::2] = phi_c
   phi[1::2,1::2] = phi_c
   return phi

def multigrid(source, bcmask):
   NITER1, NITER2, RELAXC, RELAX = 2, 5, 0.86, 0.9
   phi = zeros(source.shape, dtype=float)
   if bcmask.all():
      return phi
   # iteration going down
   for i in range(NITER1):
      res = resid(phi, source, bcmask)
      phi -= res / 4.0
   # coarse grid
   if min(source.shape) > 2:
      n = int(log2(min(source.shape))) - 2
      res = resid(phi, source, bcmask)
      res_c, bcmask_c = collcect(res, bcmask)
      phi += distribute(multigrid(res_c, bcmask_c)) * RELAXC
   phi[bcmask] = 0.0
   # iteration going up
   for i in range(NITER2):
      res = resid(phi, source, bcmask)
      phi -= res / 4.0 * RELAX
   # W-cycle, second coarse grid
   if min(source.shape) > 2:
      n = int(log2(min(source.shape))) - 2
      res = resid(phi, source, bcmask)
      res_c, bcmask_c = collcect(res, bcmask)
      phi += distribute(multigrid(res_c, bcmask_c)) * RELAXC
   phi[bcmask] = 0.0
   # W-cycle, second iteration going up
   for i in range(NITER2):
      res = resid(phi, source, bcmask)
      phi -= res / 4.0 * RELAX
   return phi

def poisson(source, bcmask, bcvalue, rel_res=1.0E-9):
   # poisson solver
   assert isinstance(source, numpy.ndarray)
   assert isinstance(bcmask, numpy.ndarray)
   assert isinstance(bcvalue, numpy.ndarray)
   assert source.ndim == bcmask.ndim == bcvalue.ndim == 2
   assert source.shape == bcmask.shape == bcvalue.shape
   assert source.dtype == float
   assert bcmask.dtype == bool
   assert bcvalue.dtype == float
   nx, ny = source.shape
   
   # initialization
   phi = zeros([nx, ny], dtype=float)

   phi[bcmask] = bcvalue[bcmask]
   res = resid(phi, source, bcmask)
   resmag0 = (res**2).sum()
   resmag = resmag0

   niter = 0
   while resmag / resmag0 > rel_res:
      res = resid(phi, source, bcmask)
      phi += multigrid(res, bcmask)
      res = resid(phi, source, bcmask)
      phi -= res / 4.0

      niter += 1
      resmag = (res**2).sum()
      print niter, resmag / resmag0

   return phi
