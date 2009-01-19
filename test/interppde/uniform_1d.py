import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from arbinterp import interp_1d
from compinterp import lagrange_1d, linear_1d, spline_1d

CASE = 'runge'
def func(x):
  f = 1.0 / (1.0 + x**2)
  return f

pylab.matplotlib.interactive(False)
n = 10
x = linspace(-5, 5, n)
y = func(x)
y = interp_1d(y, x, fx)

