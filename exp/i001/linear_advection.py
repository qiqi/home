import os

import numpy
import pylab

import rungekutta
from interp1d import interp_1d_grad_coef

# pylab.matplotlib.interactive(True)
n = 50
s = numpy.linspace(-1, 1, n+1)
dx = s[1] - s[0]
gamma = 2.0 / dx
dt = 0.005
stages = rungekutta.runge6
nstage = len(stages.alpha)

# interpolation based derivative
# x = numpy.sin(s * numpy.pi / 2)
x = s
da = numpy.zeros([n+1, n+1])
for i in range(1,n):
  if i < n - i:
    da[:2*i+1,i] = interp_1d_grad_coef(i,s[:2*i+1],gamma)
  else:
    # da[2*i-n:,i] = interp_1d_grad_coef(n-i,s[2*i-n:],gamma)
    da[:,i] = interp_1d_grad_coef(i,s,gamma)
  # da[:,i] /= numpy.cos(s[i] * numpy.pi / 2)
pylab.plot(x, da, '+-')
# pylab.show()
# stop

u = numpy.exp(-x**2*10)
u = numpy.cos(-(x+1)*10)
u0 = numpy.zeros(n+1)
du = numpy.zeros([n+1, nstage])
for i in range(1001):
  if i % 20 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-2.0,2.0])
    pylab.title('t = %f' % (i * dt))
    #pylab.savefig('advect_interp_%04d.png' % int(i*dt*10/2))
    pylab.savefig('advect_interp_%04d.png' % i)
  u0[:] = u
  du[:] = 0
  for istage in range(nstage):
    for j in range(1, n):
      du[j,istage] = -numpy.dot(da[:,j], u)
    u[:] = u0
    for jstage in range(istage + 1):
      u += du[:,jstage] * dt * stages.beta[istage][jstage]
    u[0] = numpy.cos(i*dt*10)
    u[0] = 0.0
# os.system('mencoder mf://advect_interp_*.png -mf type=png:fps=10 -ovc lavc ' + \
#           '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o advect_interp.avi')
# os.system('rm -rf advect_interp_*.png')

