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

# # second order central difference
# u = numpy.ones(n+1)
# u[x < 0] = 0.0
# u[-1] = u[0]
# u0 = numpy.zeros(n+1)
# uu = numpy.zeros(3*n+1)
# du = numpy.zeros([n+1, nstage])
# for i in range(2001):
#   if i % 20 == 0:
#     print '%d / %d' % (i, 20000)
#     pylab.cla()
#     pylab.plot(x,u,'+-')
#     pylab.ylim([-0.5,1.5])
#     pylab.title('t = %f' % (i * dt))
#     pylab.savefig('centralupwind_%04d.png' % int(i*dt*10/2))
#   u0[:] = u
#   for istage in range(nstage):
#     uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
#     du[:,istage] = -(uu[n+1:2*n+2] - uu[n-1:2*n]) / (2.0*dx)
#     u[:] = u0
#     for jstage in range(istage + 1):
#       u += du[:,jstage] * dt * stages.beta[istage][jstage]
# os.system('mencoder mf://centralupwind*.png -mf type=png:fps=10 -ovc lavc ' + \
#           '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o centralupwind.avi')
# os.system('rm -rf centralupwind_*.png')


# interpolation based derivative
xx = numpy.zeros(2*n+1)
xx[:n] = x[:-1]
xx[n:] = x + x[-1] - x[0]
gamma = 2.0 / dx
da = interp_1d_grad_coef(n,xx[:n+6],gamma)

u = numpy.ones(n+1)
u[x < 0] = 0.0
u[-1] = u[0]
u0 = numpy.zeros(n+1)
uu = numpy.zeros(3*n+1)
du = numpy.zeros([n+1, nstage])
for i in range(20001):
  if i % 20 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-0.5,1.5])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('interpupwind_%04d.png' % int(i*dt*10/2))
  u0[:] = u
  for istage in range(nstage):
    uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
    for j in range(n+1):
      du[j,istage] = -numpy.dot(da, uu[j:n+6+j])
    assert du[0,istage] == du[-1,istage]
    u[:] = u0
    for jstage in range(istage + 1):
      u += du[:,jstage] * dt * stages.beta[istage][jstage]
os.system('mencoder mf://interpupwind_*.png -mf type=png:fps=10 -ovc lavc ' + \
          '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o interpupwind.avi')
os.system('rm -rf interpupwind_*.png')

