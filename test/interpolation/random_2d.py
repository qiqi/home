import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from arbinterp import interp_nd

# CASE = 'cos_2d'
# def func(x):
#   f = numpy.cos(x.sum(1))
#   return f
# 
# CASE = 'runge_2d'
# def func(x):
#   f = 1.0 / (1.0 + (x**2).sum(1))
#   return f
# 
CASE = 'notch_2d'
def func(x):
  f = numpy.cos(x).prod(1)
  f -= 2*exp(-(x[:,1]*4)**2 - x[:,0]**2)
  return f

# CASE = 'step_2d'
# def func(x):
#   f = exp(-(x**2).sum(1))
#   f[x[:,0]<0.25*x[:,1]**2] *= -1
#   return f

def cartesian_2d(x1d):
  assert x1d.ndim == 1
  n = x1d.size
  x = zeros([n**2, 2])
  for i in range(n):
    x[i*n:i*n+n,0] = x1d
    x[i::n,1] = x1d
  return x

pylab.matplotlib.interactive(False)
n_list = numpy.arange(4, 15, 2)**2
nz = 20; z1d = linspace(-2, 2, nz); z = cartesian_2d(z1d); fze = func(z)

cc = linspace(fze.min()-0.5*(fze.max()-fze.min()), \
              fze.max()+0.5*(fze.max()-fze.min()), 30)

pylab.figure(figsize=(5,5))
pylab.contour(z1d, z1d, fze.reshape([nz,nz]), cc)
pylab.axis('scaled'); pylab.axis([-2.1, 2.1, -2.1, 2.1]);
pylab.savefig('output/%s/%s_exact.png' % (CASE, CASE))

# my interpolation
linf, l2 = [], []

xall = []
for line in file('quasi-random/niederreiter_2d.txt').readlines():
  xall.append(map(float, line.strip().split()))
xall = numpy.array(xall) * 2 - 1
# xall = 2./3 * (xall + 0.5*xall**3)
xall *= 2

for n in n_list:
  x = xall[:n,:]
  fx = func(x)
  fz, sig = interp_nd(z, x, fx, compute_dfz=True)
  pylab.figure(figsize=(5,5))
  pylab.plot(x[:,0], x[:,1], 'ok', mfc=None)
  pylab.contour(z1d, z1d, fz.reshape([nz,nz]), cc)
  pylab.axis('scaled'); pylab.axis([-2.1, 2.1, -2.1, 2.1]);
  pylab.savefig('output/%s/random_%s_1d_%03d.png' % (CASE, CASE, n))
  linf.append(abs(fz-fze).max())
  l2.append(sqrt(((fz-fze)**2).sum() / z.size))
  print n, linf[-1], l2[-1]

pylab.figure(figsize=(5,3))
pylab.semilogy(n_list, linf, ':ok')
pylab.semilogy(n_list, l2, '-ok')
pylab.savefig('output/%s/%s_convergence_random.png' % (CASE, CASE))

