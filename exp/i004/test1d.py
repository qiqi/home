"""
Qiqi Wang  January 2009  Sample use of interp1d.py
"""

import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from wanginterp import Interp1D, Interp1DVG

n = 4  # number of data points

y = linspace(-5, 5, 250)
fy = y
x = linspace(-5, 5, n)
fx = x

interp = Interp1DVG(x, fx, l=1, verbose=1)
fp, sig = interp.grad(y, compute_dfp=True)
print fp.max(), fp.min()

pylab.plot(y, fp)
pylab.plot(y, sig)
pylab.show()
