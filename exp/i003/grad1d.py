"""
Qiqi Wang  January 2009  Sample use of interp1d.py
"""

import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from wang import Interp1D, Interp1DVG

# ========================= Tune this =========================

CASE = 'cos'
def func(x):
  f = numpy.cos(x)
  return f
def grad(x):
  fp = -numpy.sin(x)
  return fp

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
#   f = numpy.cos(x)
#   f[x<0] *= -1
#   return f
# def grad(x):
#   fp = -numpy.sin(x)
#   fp[x<0] *= -1
#   return fp

n = 26  # number of data points

# ========================= Tune this =========================

y = linspace(-5, 5, 250)
fye = func(y)
fpe = grad(y)
x = linspace(-5, 5, n)
fx = func(x)
fp = grad(x)

# --------- use of Interp1D class ----------
interp = Interp1D(x, fx, l=1, verbose=1)
fy, sig = interp.grad(y, compute_dfp=True)
# ------------------------------------------

pylab.figure(figsize=(10,6))
pylab.plot(y, fpe, '--r')
pylab.plot(y, fy, 'k')
pylab.plot(y, fy + 3*sig, ':b')
pylab.plot(y, fy - 3*sig, ':b')
pylab.plot(x, fp, 'ok')
pylab.ylim([fpe.min()-0.2*(fpe.max()-fpe.min()), \
            fpe.max()+0.2*(fpe.max()-fpe.min())])
pylab.xlabel('$x$', verticalalignment='center', fontsize=16);
pylab.ylabel('$f(x)$', horizontalalignment='right', fontsize=16)
pylab.savefig('1g.png')

linf = abs(fy-fpe).max()
l2 = sqrt(((fy-fpe)**2).sum() / y.size)

print
print "red dash line: target function gradient"
print "black solid line: numerical gradient"
print "dotted blue lines: 3 sigma confidence interval"
print
print "Number of points: %d;  L-infinity error: %f;  L-2 error: %f" % \
      (n, linf, l2)
print

# --------- use of Interp1DVG class ----------
interp = Interp1DVG(x, fx, l=1, verbose=1)
fy, sig = interp.grad(y, compute_dfp=True)
# ------------------------------------------

pylab.figure(figsize=(10,6))
pylab.plot(y, fpe, '--r')
pylab.plot(y, fy, 'k')
pylab.plot(y, fy + 3*sig, ':b')
pylab.plot(y, fy - 3*sig, ':b')
pylab.plot(x, fp, 'ok')
pylab.ylim([fpe.min()-0.2*(fpe.max()-fpe.min()), \
            fpe.max()+0.2*(fpe.max()-fpe.min())])
pylab.xlabel('$x$', verticalalignment='center', fontsize=16);
pylab.ylabel('$f(x)$', horizontalalignment='right', fontsize=16)
pylab.savefig('2g.png')

linf = abs(fy-fpe).max()
l2 = sqrt(((fy-fpe)**2).sum() / y.size)

print
print "red dash line: target function gradient"
print "black solid line: numerical gradient"
print "dotted blue lines: 3 sigma confidence interval"
print
print "Number of points: %d;  L-infinity error: %f;  L-2 error: %f" % \
      (n, linf, l2)
print

