import numpy
import pylab
from numpy import linspace, exp, sqrt, random, pi

from arbinterp import interp_1d_adaptive
from compinterp import floater_hormann_adaptive

pylab.matplotlib.interactive(False)

y = linspace(0, 1, 250)
fye = exp(-((y-0.5)*16)**2)

n_list = range(4, 50, 4)
randx = random.random(size=n_list[-1])
randx = numpy.tan((randx-0.5)*pi*0.6) / numpy.tan(0.3*pi) * 0.5 + 0.5

linf, l2 = [], []
linf_m, l2_m = [], []
for n in n_list:
  # x = linspace(0, 1, n)
  x = numpy.array(sorted(randx[:n]))
  fx = exp(-((x-0.5)*16)**2)
  fy = floater_hormann_adaptive(y, x, fx)
  fy_m = interp_1d_adaptive(y, x, fx)
  pylab.figure()
  pylab.plot(y, fye, '--r')
  pylab.plot(y, fy, 'g')
  pylab.plot(y, fy_m, 'b')
  pylab.plot(x, fx, '+k')
  pylab.ylim([-0.2, 1.2])
  linf.append(abs(fy-fye).max())
  l2.append(sqrt(((fy-fye)**2).sum()))
  linf_m.append(abs(fy_m-fye).max())
  l2_m.append(sqrt(((fy_m-fye)**2).sum()))
  print linf[-1], l2[-1], ' mine ', linf_m[-1], l2_m[-1]

pylab.figure()
pylab.semilogy(n_list, linf, '-+g')
pylab.semilogy(n_list, l2, '-+g')
pylab.semilogy(n_list, linf_m, '-+b')
pylab.semilogy(n_list, l2_m, '-+b')

pylab.show()
