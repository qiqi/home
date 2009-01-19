import os

import numpy
import pylab

import rungekutta
from arbinterp import interp_1d_grad_coef

n = 50
x = numpy.linspace(-1, 1, n+1)
dx = x[1] - x[0]
dt = 0.01
stages = rungekutta.runge6
nstage = len(stages.alpha)

# second order central difference
u = numpy.exp(-25*x**2)
u0 = numpy.zeros(n+1)
uu = numpy.zeros(3*n+1)
du = numpy.zeros([n+1, nstage])
for i in range(4001):
  if i % 20 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-0.5,1.0])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('central_%04d.png' % int(i*dt*10/2))
  u0[:] = u
  for istage in range(nstage):
    uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
    du[:,istage] = -(uu[n+1:2*n+2] - uu[n-1:2*n]) / (2.0*dx)
    u[:] = u0
    for jstage in range(istage + 1):
      u += du[:,jstage] * dt * stages.beta[istage][jstage]
os.system('mencoder mf://central_*.png -mf type=png:fps=10 -ovc lavc ' + \
          '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o central.avi')
os.system('rm -rf central_*.png')


# interpolation based derivative
gamma = 2.0 / dx
xx = numpy.zeros(2*n+1)
xx[:n] = x[:-1]
xx[n:] = x + x[-1] - x[0]
da = interp_1d_grad_coef(n,xx,gamma)

u = numpy.exp(-25*x**2)
u0 = numpy.zeros(n+1)
uu = numpy.zeros(3*n+1)
du = numpy.zeros([n+1, nstage])
for i in range(20001):
  if i % 20 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-0.5,1.0])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('interp_%04d.png' % int(i*dt*10/2))
  u0[:] = u
  for istage in range(nstage):
    uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
    for j in range(n+1):
      du[j,istage] = -numpy.dot(da, uu[j:2*n+1+j])
    u[:] = u0
    for jstage in range(istage + 1):
      u += du[:,jstage] * dt * stages.beta[istage][jstage]
os.system('mencoder mf://interp_*.png -mf type=png:fps=10 -ovc lavc ' + \
          '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o interp.avi')
os.system('rm -rf interp_*.png')

