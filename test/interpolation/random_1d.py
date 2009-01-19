import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from interp1d import interp_1d
from compinterp import lagrange_1d, linear_1d, spline_1d

CASE = 'cos'
def func(x):
  f = numpy.cos(x)
  return f

# CASE = 'runge'
# def func(x):
#   f = 1.0 / (1.0 + x**2)
#   return f
# 
# CASE = 'notch'
# def func(x):
#   f = numpy.cos(x)
#   f -= 2*exp(-(x*4)**2)
#   return f
# 
# CASE = 'step'
# def func(x):
#   f = exp(-x**2)
#   f[x<0] *= -1
#   return f

pylab.matplotlib.interactive(False)
n_list = range(8, 128, 8)
y = linspace(-5, 5, 250)
fye = func(y)

pylab.figure(figsize=(5,3))
pylab.plot(y, fye, '-k')
pylab.ylim([fye.min()-0.2*(fye.max()-fye.min()), \
            fye.max()+0.2*(fye.max()-fye.min())])
pylab.savefig('output/%s/%s_exact.png' % (CASE, CASE))

# my interpolation
linf, l2 = [], []

xall = []
for line in file('quasi-random/niederreiter_1d.txt').readlines():
  xall.append(float(line))
xall = numpy.array(xall) * 10 - 5

for n in n_list:
  x = xall[:n]
  fx = func(x)
  fy, sig = interp_1d(y, x, fx, compute_dfz=True)
  pylab.figure(figsize=(5,3))
  pylab.plot(y, fye, '--k')
  pylab.plot(y, fy, 'k')
  # pylab.plot(y, fy + sig, ':b')
  # pylab.plot(y, fy - sig, ':b')
  # pylab.plot(y, fy_l, '-.g')
  pylab.plot(x, fx, 'ok')
  margin = 0.2
  if CASE == 'step': margin = 0.5
  pylab.ylim([fye.min()-margin*(fye.max()-fye.min()), \
              fye.max()+margin*(fye.max()-fye.min())])
  pylab.xlabel('$x$', verticalalignment='center', fontsize=16);
  pylab.ylabel('$f(x)$', horizontalalignment='right', fontsize=16)
  pylab.savefig('output/%s/random_%s_1d_%03d.png' % (CASE, CASE, n))
  linf.append(abs(fy-fye).max())
  l2.append(sqrt(((fy-fye)**2).sum() / y.size))
  print n, linf[-1], l2[-1]

pylab.figure(figsize=(5,3))
pylab.semilogy(n_list, linf, ':ok')
pylab.semilogy(n_list, l2, '-ok')
pylab.ylabel('approximation error', horizontalalignment='center', fontsize=16)
pylab.savefig('output/%s/%s_convergence_random.png' % (CASE, CASE))

# spline
linf_s, l2_s = [], []
for n in n_list:
  x = xall[:n]
  fx = func(x)
  fy = spline_1d(y, x, fx)
  linf_s.append(abs(fy-fye).max())
  l2_s.append(sqrt(((fy-fye)**2).sum() / y.size))

pylab.figure(figsize=(5,3))
pylab.semilogy(n_list, linf, '-ok')
pylab.semilogy(n_list, l2, ':ok')
pylab.semilogy(n_list, linf_s, '-+k', markersize=8)
pylab.semilogy(n_list, l2_s, ':+k', markersize=8)
pylab.savefig('output/%s/%s_convergence_random_spline.png' % (CASE, CASE))

# pylab.show()

