import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, dot, exp, sqrt, cos, pi

from arbinterp import interp_1d, interp_1d_adaptive
from compinterp import lagrange_1d

pylab.matplotlib.interactive(False)
gamma = 16.0
y = linspace(0, 1, 250)
fye = exp(-((y-0.5)*16)**2)

n_list = range(4, 50, 2)
linf, l2 = [], []
linf_l, l2_l = [], []

for n in n_list:
  x = (1 - cos(linspace(0, pi, n))) / 2.0
  fx = exp(-((x-0.5)*16)**2)
  # fy, sig = interp_1d(y, x, fx, gamma, compute_sigma=True)
  fy, sig = interp_1d_adaptive(y, x, fx, compute_sigma=True)
  fy_l = lagrange_1d(y, x, fx)
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
  linf_l.append(abs(fy_l-fye).max())
  l2_l.append(sqrt(((fy_l-fye)**2).sum()))
  print n, linf[-1], l2[-1], linf_l[-1], l2_l[-1]


pylab.figure()
pylab.semilogy(n_list, linf, '-+b')
pylab.semilogy(n_list, l2, '-+b')
pylab.semilogy(n_list, linf_l, '-+g')
pylab.semilogy(n_list, l2_l, '-+g')
pylab.show()

