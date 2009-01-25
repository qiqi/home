"""
Qiqi Wang  January 2009  Sample use of interp1d.py
"""

import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from wanginterp import Interp2D

# def func(x):
#   return x[:,0]**2 + x[:,1]**2
# 
# def func(x):
#   return x[:,0] + x[:,1]
# 
def func(x):
  return 1.0 / (x[:,0]**2 + x[:,1]**2 + 1.0)

xall = []
for line in file('quasi-random/niederreiter_2d.txt').readlines():
  xall.append(map(float, line.strip().split()))
xall = numpy.array(xall) * 2 - 1
# xall = 2./3 * (xall + 0.5*xall**3)
xall *= 2

n = 36  # number of data points
x = xall[:n,:]


def cartesian_2d(x1d):
  assert x1d.ndim == 1
  n = x1d.size
  x = zeros([n**2, 2])
  for i in range(n):
    x[i*n:i*n+n,0] = x1d
    x[i::n,1] = x1d
  return x

ny = 20
y1d = linspace(-2, 2, ny)
y = cartesian_2d(y1d)

fy = func(y)
fx = func(x)

interp = Interp2D(x, fx, l=1, verbose=1)
f, sig = interp.interp(y, compute_df=True)

pylab.figure()
c1 = pylab.contourf(y1d, y1d, f.reshape([ny,ny]))
pylab.colorbar(c1, fraction=0.02)
pylab.plot(x[:,0], x[:,1], 'ok')
pylab.axis('scaled')

pylab.figure()
c1 = pylab.contourf(y1d, y1d, sig.reshape([ny,ny]))
pylab.colorbar(c1, fraction=0.02)
pylab.plot(x[:,0], x[:,1], 'ok')
pylab.axis('scaled')

pylab.figure()
c1 = pylab.contourf(y1d, y1d, abs(f-fy).reshape([ny,ny]))
pylab.colorbar(c1, fraction=0.02)
pylab.plot(x[:,0], x[:,1], 'ok')
pylab.axis('scaled')

pylab.show()

