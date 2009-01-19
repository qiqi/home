from arbinterp import *
from numpy import *
import pylab

xi = -5 + 10./4
dye = -2.0 * xi / (1 + xi**2)**2
ns, es, es2 = [], [], []
for n in range(5,80,4):
  x = linspace(-5,5,n)
  y = 1 / (1 + x**2)
  z = x[1]
  # gamma = calc_gamma_1d(x, y, zeros(n), zeros(0), zeros(0), zeros(0))
  gamma = 1.8/(x[1] - x[0])
  ni = int((n-1)/4)
  assert x[ni] == xi
  da = -interp_1d_grad_coef(ni,x,gamma)
  dy = dot(y, da)
  ns.append(n)
  es.append(abs(dy - dye))
  es2.append(abs((y[ni+1] - y[ni-1]) / (x[ni+1] - x[ni-1]) - dye))

pylab.semilogy(ns, es, '+-')
pylab.semilogy(ns, es2, 'x:')
pylab.show()

# n = 45
# ix = 20
# gamma = 3
# x = linspace(-5,5,n)
# da = interp_1d_grad_coef(ix,x,gamma)
# nil = zeros(0)
# a0, b, er2 = interp_1d_coef(x[ix], x, zeros(n), nil, nil, 1, gamma)
# a1, b, er2 = interp_1d_coef(x[ix] + 0.01, x, zeros(n), nil, nil, 1, gamma)
# pylab.plot(-da)
# pylab.plot((a1 - a0) / 0.01, ':')
# pylab.show()
