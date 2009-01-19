import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from interp1d import interp_1d_coef

pylab.matplotlib.interactive(False)

gamma = [1, 10, 25, 1000]
z = linspace(-0.2, 1.2, 250)
l = 1

for n in [7]:
  x = linspace(0, 1, n)
  pylab.figure(figsize=(10,5))
  for i in range(len(gamma)):
    a = zeros([z.size,n]); er2 = zeros(z.size)
    for j in range(z.size):
      a[j,:], b, er2[j] = interp_1d_coef(z[j], x, zeros(n), zeros(0), zeros(0),\
                                         1.0, gamma[i], n, l)
    pylab.subplot(2,2,i+1)
    for j in range(n):
      pylab.plot(z, a[:,j])
    for j in range(n):
      pylab.plot([x[j], x[j]], [-0.8,1.8], ':')
    pylab.xlim([z.min(), z.max()])
    pylab.ylim([-0.8,1.8])
    pylab.xlabel('$x$', verticalalignment='center', fontsize=16);
    pylab.ylabel('$a_i(x)$', horizontalalignment='center', fontsize=16);

pylab.savefig('basis_uniform.png')

for n in [7]:
  x = zeros(n)
  while (x[1:] - x[:-1]).min() < 0.05:
    x = numpy.array([0] + sorted(numpy.random.random(size=n-2)) + [1])
  pylab.figure(figsize=(10,5))
  for i in range(len(gamma)):
    a = zeros([z.size,n])
    for j in range(z.size):
      a[j,:], b, er2 = interp_1d_coef(z[j], x, zeros(n), zeros(0), zeros(0), \
                                      1.0, gamma[i], n, l)
    pylab.subplot(2,2,i+1)
    for j in range(n):
      pylab.plot(z, a[:,j])
    for j in range(n):
      pylab.plot([x[j], x[j]], [-0.8,1.8], ':')
    pylab.xlim([-0.2,1.2])
    pylab.ylim([-0.8,1.8])
    pylab.xlabel('$x$', verticalalignment='center', fontsize=16);
    pylab.ylabel('$a_i(x)$', horizontalalignment='center', fontsize=16);

pylab.savefig('basis_random.png')
pylab.show()
