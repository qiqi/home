import numpy
import pylab
from numpy import zeros, logical_or, arange, absolute, log

from laplace import poisson

resids, err = [], []
NGRID = [8, 16, 32, 64, 128, 256, 512, 1024]
for N in NGRID:
   source = zeros([N,N], dtype=float)
   x = zeros([N,N], dtype=float)
   y = zeros([N,N], dtype=float)
   x[:] = (arange(N) / float(N-1) - 0.5) * 8.0
   y.transpose()[:] = (arange(N) / float(N-1) - 0.5) * 8.0
   bcmask = logical_or(x**2 + y**2 >= 16.0, x**2 + y**2 <= 1.0)
   bcvalue = zeros([N,N], dtype=float)
   bcvalue[x**2 + y**2 <= 1.0] = 1.0
   
   phi, relres = poisson(source, bcmask, bcvalue)
   pylab.figure(); pylab.contour(x, y, phi)
   pylab.axis('scaled'); pylab.savefig('soln%04d.png' % N)
   resids.append(relres)

   phi_exact = (log(16.0) - log(x**2 + y**2)) / log(16.0)
   phi_exact[bcmask] = phi[bcmask]
   err.append(absolute(phi - phi_exact).max())

pylab.figure()
legends = []
for relres, N in zip(resids, NGRID):
   pylab.semilogy(range(0, len(relres)*2, 2), relres, '+-')
   legends.append('%d by %d' % (N, N))
pylab.semilogy([0, 18], [1.0E-9, 1.0E-9], ':')
pylab.legend(legends)
pylab.savefig('conv.png')

pylab.figure()
pylab.loglog(NGRID, err, '+-')
pylab.savefig('order.png')
