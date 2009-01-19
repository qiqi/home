import numpy
import pylab
from numpy import ones, zeros, linspace, kron, dot, exp, sqrt, diag, pi, tan

from arbinterp import interp_2d

def lattice2d(n):
  xy = linspace(0, 1, n)
  xy = kron(xy.reshape(n,1), ones([1,n]))
  y = numpy.array([xy.reshape(n**2), xy.transpose().reshape(n**2)])
  return y.transpose()

gamma = 16.0
ny = 10; y = lattice2d(ny) * 0.2 + 0.5
yout = lattice2d(ny) * 0.3 + 0.8
y_half = numpy.array(filter(lambda a: a[0]>=a[1], y))
yout_half = numpy.array(filter(lambda a: a[0]>=a[1], yout))
y0 = numpy.array([0.5, 0.5])
fye = exp(-(((y-y0)*16)**2).sum(1))

n_list = range(5, 18, 2)
linf, l2 = [], []

for n in n_list:
  x = lattice2d(n)
  x = tan((x-0.5)*pi*0.8) / tan(0.4*pi) * 0.5 + 0.5
  # x = (1 - cos(x * pi)) / 2.0
  fx = exp(-(((x-y0)*16)**2).sum(1))
  fy_half = interp_2d(y_half, x, fx, gamma)
  fyout_half = interp_2d(yout_half, x, fx, gamma)
  fy, fyout = zeros(y.shape[0]), zeros(y.shape[0])
  for i in range(len(fy)):
    ix, iy = max(i/ny, i%ny), min(i/ny, i%ny)
    fy[i] = fy_half[ix*(ix+1)/2 + iy]
    fyout[i] = fyout_half[ix*(ix+1)/2 + iy]
  pylab.figure()
  z = linspace(-0.1,1.1,10)
  pylab.contour(y[:ny,1], y[:ny,1], fy.reshape(ny,ny), z)
  pylab.contour(1-y[:ny,1], y[:ny,1], fy.reshape(ny,ny), z)
  pylab.contour(1-y[:ny,1], 1-y[:ny,1], fye.reshape(ny,ny), z)
  pylab.contour(y[:ny,1], 1-y[:ny,1], fye.reshape(ny,ny), z)
  pylab.plot(x[:,0], x[:,1], '+k')
  pylab.axis('scaled')
  pylab.axis([-0.1, 1.1, -0.1, 1.1])
  linf.append(abs(fy-fye).max())
  l2.append(sqrt(((fy-fye)**2).sum()))
  print n, linf[-1], l2[-1], abs(fyout).max()

pylab.figure()
pylab.semilogy(n_list, linf, '-+')
pylab.semilogy(n_list, l2, '-+')
pylab.show()
