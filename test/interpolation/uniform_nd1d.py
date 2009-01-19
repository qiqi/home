import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from arbinterp import interp_nd

# CASE = 'cos'
# def func(x):
#   f = numpy.cos(x)
#   return f
# 
# CASE = 'runge'
# def func(x):
#   f = 1.0 / (1.0 + x**2)
#   return f
# 
CASE = 'notch'
def func(x):
  f = numpy.cos(x)
  f -= 2*exp(-(x*4)**2)
  return f

# CASE = 'step'
# def func(x):
#   f = exp(-x**2)
#   f[x<0] *= -1
#   return f

pylab.matplotlib.interactive(False)
n_list = range(8, 128, 8)
y = linspace(-5, 5, 250).reshape([250,1])
fye = func(y).reshape([250])

pylab.figure(figsize=(5,3))
pylab.plot(y, fye, '-k')
pylab.ylim([fye.min()-0.2*(fye.max()-fye.min()), \
            fye.max()+0.2*(fye.max()-fye.min())])
pylab.savefig('output/%s/%s_exact.png' % (CASE, CASE))

# my interpolation
linf, l2 = [], []
for n in n_list:
  x = linspace(-5, 5, n).reshape([n,1])
  fx = func(x).reshape([n])
  fy, sig = interp_nd(y, x, fx, compute_dfz=True)
  fy = fy.reshape([fy.size])
  pylab.figure(figsize=(5,3))
  pylab.plot(y, fye, '--k')
  pylab.plot(y, fy, 'k')
  # pylab.plot(y, fy + sig, ':b')
  # pylab.plot(y, fy - sig, ':b')
  # pylab.plot(y, fy_l, '-.g')
  pylab.plot(x, fx, 'ok')
  pylab.ylim([fye.min()-0.2*(fye.max()-fye.min()), \
              fye.max()+0.2*(fye.max()-fye.min())])
  pylab.savefig('output/%s/uniform_%s_nd1d_%03d.png' % (CASE, CASE, n))
  linf.append(abs(fy-fye).max())
  l2.append(sqrt(((fy-fye)**2).sum() / y.size))
  print n, linf[-1], l2[-1]

pylab.figure(figsize=(5,3))
pylab.semilogy(n_list, linf, ':ok')
pylab.semilogy(n_list, l2, '-ok')
pylab.savefig('output/%s/%s_convergence_uniform_nd1d.png' % (CASE, CASE))

