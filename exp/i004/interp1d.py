"""
Qiqi Wang  January 2009  Sample use of interp1d.py
"""

import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from wanginterp import Interp1D, Interp1DVG

# ========================= Tune this =========================

#CASE = 'cos'
#def func(x):
#  f = numpy.cos(x)
#  return f
#
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
CASE = 'step'
def func(x):
  f = numpy.cos(x)
  f[x<0] *= -1
  return f

n = 16  # number of data points

# ========================= Tune this =========================

y = linspace(-5, 5, 250)
fye = func(y)
x = linspace(-5, 5, n)
fx = func(x)

# --------- use of Interp1D class ----------
interp = Interp1D(x, fx, l=1, verbose=1)
fy, sig = interp.interp(y, compute_df=True)
# ------------------------------------------

pylab.figure(figsize=(10,6))
pylab.plot(y, fye, '--r')
pylab.plot(y, fy, 'k')
pylab.plot(y, fy + 3*sig, ':b')
pylab.plot(y, fy - 3*sig, ':b')
pylab.plot(x, fx, 'ok')
pylab.ylim([fye.min()-0.2*(fye.max()-fye.min()), \
            fye.max()+0.2*(fye.max()-fye.min())])
pylab.xlabel('$x$', verticalalignment='center', fontsize=16);
pylab.ylabel('$f(x)$', horizontalalignment='right', fontsize=16)
pylab.savefig('1.png')

linf = abs(fy-fye).max()
l2 = sqrt(((fy-fye)**2).sum() / y.size)

print
print "red dash line: target function"
print "black solid line: interpolant"
print "dotted blue lines: 3 sigma confidence interval"
print
print "Number of points: %d;  L-infinity error: %f;  L-2 error: %f" % \
      (n, linf, l2)
print

# --------- use of Interp1DVG class ----------
interp = Interp1DVG(x, fx, l=1, verbose=1)
fy, sig = interp.interp(y, compute_df=True)
# ------------------------------------------

pylab.figure(figsize=(10,6))
pylab.plot(y, fye, '--r')
pylab.plot(y, fy, 'k')
pylab.plot(y, fy + 3*sig, ':b')
pylab.plot(y, fy - 3*sig, ':b')
pylab.plot(x, fx, 'ok')
pylab.ylim([fye.min()-0.2*(fye.max()-fye.min()), \
            fye.max()+0.2*(fye.max()-fye.min())])
pylab.xlabel('$x$', verticalalignment='center', fontsize=16);
pylab.ylabel('$f(x)$', horizontalalignment='right', fontsize=16)
pylab.savefig('2.png')

linf = abs(fy-fye).max()
l2 = sqrt(((fy-fye)**2).sum() / y.size)

print
print "red dash line: target function"
print "black solid line: interpolant"
print "dotted blue lines: 3 sigma confidence interval"
print
print "Number of points: %d;  L-infinity error: %f;  L-2 error: %f" % \
      (n, linf, l2)
print

