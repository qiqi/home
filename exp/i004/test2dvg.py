"""
Qiqi Wang  January 2009  Sample use of interp1d.py
"""

import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi, log10, asarray, floor

from wanginterp import Interp2DVG

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

n = 56  # number of data points
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

interp = Interp2DVG(x, fx, l=1, verbose=1, safety_factor=2.0)
f, sig = interp.interp(y, compute_df=True)
gamma = exp(interp.log_gamma_interp.interp(y))

def contourplot(y1d, f, title):
  mag = 10**int(floor(log10((f.max() - f.min()) / 10.)))
  fmin, fmax = f.min() - mag, f.max() + mag
  if f.min() >= 0:
    fmin = max(fmin, 0)
  levels = mag * asarray(linspace(fmin, fmax, 10) / mag, int)
  c3 = pylab.contourf(y1d, y1d, f.reshape([len(y1d), len(y1d)]), levels)
  b3 = pylab.colorbar(c3, fraction=0.02)
  pylab.plot(x[:,0], x[:,1], 'ow')
  pylab.title(title)
  pylab.axis('scaled')

pylab.figure(figsize=(8,6))

pylab.subplot(2,2,1)
contourplot(y1d, f, 'interpolant')

pylab.subplot(2,2,2)
contourplot(y1d, sig, 'estimated residual')

pylab.subplot(2,2,3)
contourplot(y1d, abs(f-fy), 'true residual')

pylab.subplot(2,2,4)
contourplot(y1d, gamma, 'gamma')

pylab.show()

