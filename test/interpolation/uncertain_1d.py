import numpy
import pylab
from numpy import zeros, ones, linspace, exp, sqrt, diag, random, \
                  arctan, tan, pi

from arbinterp import interp_1d, lagrange_1d, interp_1d_adaptive

pylab.matplotlib.interactive(False)
gamma = 36.0
y = linspace(0, 1, 100)
fye = exp(-((y-0.5)*16)**2)
# fye[:] = 0; fye[y>0.5] = 1
# fye = arctan(100*(y-0.5)) / 3.0 + 0.5

n_list = range(4, 50, 4)
linf, l2 = [], []

xall = random.random(size=n_list[-1])
dall = exp(-3.0 * random.random(size=n_list[-1])) * 0.05
eall = random.normal(size=n_list[-1]) * dall

for n in n_list:
  x = xall[:n]
  x = tan((x-0.5)*pi*0.8) / tan(0.4*pi) * 0.5 + 0.5
  fx = exp(-((x-0.5)*16)**2) + eall[:n]
  dfx = dall[:n]
  # fx[:] = 0; fx[x>0.5] = 1
  # fx = arctan(100*(x-0.5)) / 3.0 + 0.5
  fy, sig = interp_1d_adaptive(y, x, fx, dfx, compute_sigma=True)
  fy_l = lagrange_1d(y, x, fx)
  if abs(fy_l).max() > 1.0E3: fy_l[:] = 0.0
  sig[sig>1000] = 1000; sig[sig<-1000] = -1000
  pylab.figure()
  pylab.plot(y, fy, 'b')
  pylab.plot(y, fy + sig, ':b')
  pylab.plot(y, fy - sig, ':b')
  pylab.plot(y, fy_l, '-.g')
  pylab.plot(y, fye, '--r')
  pylab.bar(x, x*0, bottom=fx, yerr=dfx, width=0)
  pylab.ylim([-0.2, 1.2])
  linf.append(abs(fy-fye).max())
  l2.append(sqrt(((fy-fye)**2).sum()))
  print n, linf[-1], l2[-1]

pylab.figure()
pylab.semilogy(n_list, linf, '-+')
pylab.semilogy(n_list, l2, '-+')
pylab.show()

