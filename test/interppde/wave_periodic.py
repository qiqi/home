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
u = numpy .zeros([n+1,2])
u[:,0] = numpy.exp(-25*x**2)
u[:,1] = u[:,0]
u0 = numpy.zeros([n+1,2])
uu = numpy.zeros([3*n+1,2])
du = numpy.zeros([n+1, 2, nstage])
for i in range(1001):
  if i % 20 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-1.0,1.0])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('periodic_central_%04d.png' % int(i*dt*10/2))
  u0[:,:] = u
  for istage in range(nstage):
    uu[:n,:] = u[:-1,:]; uu[n:2*n,:] = u[:-1,:]; uu[2*n:,:] = u
    du[:,0,istage] = (uu[n+1:2*n+2,1] - uu[n-1:2*n,1]) / (2.0*dx)
    du[:,1,istage] = (uu[n+1:2*n+2,0] - uu[n-1:2*n,0]) / (2.0*dx)
    u[:,:] = u0
    for jstage in range(istage + 1):
      u += du[:,:,jstage] * dt * stages.beta[istage][jstage]
os.system('mencoder mf://periodic_central_*.png -mf type=png:fps=10 -ovc lavc '\
    + '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o periodic_central.avi')
os.system('rm -rf periodic_central_*.png')


# # # interpolation based derivative
# gamma = 2.0 / dx
# da = numpy.zeros([n+1, n+1])
# for i in range(n+1):
#   da[:,i] = interp_1d_grad_coef(i,x,gamma)
# 
# u = numpy .zeros([n+1, 2])
# u[:,0] = numpy.exp(-25*x**2)
# u[:,1] = u[:,0]
# u0 = numpy.zeros([n+1, 2])
# du = numpy.zeros([n+1, 2, nstage])
# for i in range(1001):
#   if i % 20 == 0:
#     print '%d / %d' % (i, 20000)
#     pylab.cla()
#     pylab.plot(x,u,'+-')
#     pylab.ylim([-0.5,1.0])
#     pylab.title('t = %f' % (i * dt))
#     pylab.savefig('wave_interp_%04d.png' % int(i*dt*10/2))
#   u0[:,:] = u
#   for istage in range(nstage):
#     for j in range(n+1):
#       du[j,0,istage] = numpy.dot(da[:,j], u[:,1])
#       du[j,1,istage] = numpy.dot(da[:,j], u[:,0])
#     u[:,:] = u0
#     for jstage in range(istage + 1):
#       u += du[:,:,jstage] * dt * stages.beta[istage][jstage]
#     u[0,0] = 0.0; u[-1,0] = 0.0
# os.system('mencoder mf://wave_interp_*.png -mf type=png:fps=10 -ovc lavc ' + \
#           '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o wave_interp.avi')
# os.system('rm -rf wave_interp_*.png')

