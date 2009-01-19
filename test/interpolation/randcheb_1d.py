import numpy
import pylab
from numpy import zeros, ones, linspace, dot, exp, sqrt, diag, random, cos, pi

from arbinterp import interp_1d, interp_1d_adaptive
from compinterp import lagrange_1d

pylab.matplotlib.interactive(False)
gamma = 16.0
y = linspace(0, 1, 250)
fye = exp(-((y-0.5)*16)**2)

n_list = range(4, 50, 4)
linf, l2 = [], []

xall = (cos(random.random(size=n_list[-1]) * pi) + 1) / 2.0

for n in n_list:
  x = xall[:n]
  fx = exp(-((x-0.5)*16)**2)
  fy, sig = interp_1d(y, x, fx, gamma, compute_sigma=True)
  fy_l = lagrange_1d(y, x, fx)
  if abs(fy_l).max() > 1.0E3: fy_l[:] = 0.0
  sig[sig>1000] = 1000; sig[sig<-1000] = -1000
  pylab.figure()
  pylab.plot(y, fy, 'b')
  pylab.plot(y, fy + sig, ':b')
  pylab.plot(y, fy - sig, ':b')
  pylab.plot(y, fy_l, '-.g')
  pylab.plot(y, fye, '--r')
  pylab.plot(x, fx, '+')
  pylab.ylim([-0.2, 1.2])
  linf.append(abs(fy-fye).max())
  l2.append(sqrt(((fy-fye)**2).sum()))
  print n, linf[-1], l2[-1]

pylab.figure()
pylab.semilogy(n_list, linf, '-+')
pylab.semilogy(n_list, l2, '-+')
pylab.show()

