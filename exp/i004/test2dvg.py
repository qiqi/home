"""
Qiqi Wang  January 2009  Sample use of interp1d.py
"""

import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi, log10, asarray, floor

import misc
from wanginterp import Interp2DVG

# def func(x):
#   return x[:,0]**2 + x[:,1]**2
# 
# def func(x):
#   return x[:,0] + x[:,1]
# 
def func(x):
  return 1.0 / (x[:,0]**2 + x[:,1]**2 + 1.0)

n = 16  # number of data points
x = (misc.niederreiter_2d[:n,:] * 2 - 1) * 2


ny = 20
y1d = linspace(-2, 2, ny)
y = misc.cartesian_2d(y1d)

fx, fy = func(x), func(y)

interp = Interp2DVG(x, fx, l=1, verbose=1, safety_factor=2.0)
f, sig = interp.interp(y, compute_df=True)
gamma = exp(interp.log_gamma_interp.interp(y))

pylab.figure(figsize=(8,6))
pylab.subplot(2,2,1)
misc.contourplot(y1d, f, x, 'interpolant')
pylab.subplot(2,2,2)
misc.contourplot(y1d, sig, x, 'estimated residual')
pylab.subplot(2,2,3)
misc.contourplot(y1d, abs(f-fy), x, 'true residual')
pylab.subplot(2,2,4)
misc.contourplot(y1d, gamma, x, 'gamma')

pylab.show()

