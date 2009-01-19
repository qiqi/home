import os

import numpy
import pylab

import rungekutta
from arbinterp import interp_1d_grad_coef

n = 100
x = numpy.linspace(-1, 1, n+1)
dx = x[1] - x[0]
dt = 0.01
stages = rungekutta.runge6
nstage = len(stages.alpha)

# second order central difference
u = numpy.ones(n+1)
u[x < 0] = 0.0
u[-1] = u[0]
u0 = numpy.zeros(n+1)
uu = numpy.zeros(3*n+1)
du = numpy.zeros([n+1, nstage])
for i in range(2001):
  if i % 20 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-0.5,1.5])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('centraldiff_%04d.png' % int(i*dt*10/2))
  u0[:] = u
  for istage in range(nstage):
    uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
    du[:,istage] = -(uu[n+1:2*n+2] - uu[n-1:2*n]) / (2.0*dx)
    dudu = uu[1:] - uu[:-1]
    ddu = dudu[1:] - dudu[:-1]
    dddu = ddu[1:] - ddu[:-1]
    ddddu = dddu[1:] - dddu[:-1]
    du[:,istage] += 0.002*dddu[n-2:2*n-1] / dx
    u[:] = u0
    for jstage in range(istage + 1):
      u += du[:,jstage] * dt * stages.beta[istage][jstage]
os.system('mencoder mf://centraldiff_*.png -mf type=png:fps=10 -ovc lavc ' + \
          '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o centraldiff.avi')
os.system('rm -rf centraldiff_*.png')


# interpolation based derivative
gamma = 2.0 / dx
xx = numpy.zeros(2*n+1)
xx[:n] = x[:-1]
xx[n:] = x + x[-1] - x[0]
da = interp_1d_grad_coef(n,xx,gamma)

u = numpy.ones(n+1)
u[x < 0] = 0.0
u[-1] = u[0]
u0 = numpy.zeros(n+1)
uu = numpy.zeros(3*n+1)
du = numpy.zeros([n+1, nstage])
for i in range(2001):
  if i % 20 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-0.5,1.5])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('interpdiff_%04d.png' % int(i*dt*10/2))
  u0[:] = u
  for istage in range(nstage):
    uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
    for j in range(n+1):
      du[j,istage] = -numpy.dot(da, uu[j:2*n+1+j])
    assert du[0,istage] == du[-1,istage]
    dudu = uu[1:] - uu[:-1]
    ddu = dudu[1:] - dudu[:-1]
    dddu = ddu[1:] - ddu[:-1]
    ddddu = dddu[1:] - dddu[:-1]
    dddddu = ddddu[1:] - ddddu[:-1]
    ddddddu = dddddu[1:] - dddddu[:-1]
    du[:,istage] += 0.002*dddddu[n-2:2*n-1] / dx
    u[:] = u0
    for jstage in range(istage + 1):
      u += du[:,jstage] * dt * stages.beta[istage][jstage]
os.system('mencoder mf://interpdiff_*.png -mf type=png:fps=10 -ovc lavc ' + \
          '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o interpdiff.avi')
os.system('rm -rf interpdiff_*.png')

