import os

import numpy
import pylab

import rungekutta
from interp1d import interp_1d_grad_coef

# pylab.matplotlib.interactive(True)

n = 50
x = numpy.linspace(-1, 1, n+1)
dx = x[1] - x[0]
dt = 0.01
stages = rungekutta.runge6
nstage = len(stages.alpha)

# # second order central difference
# u = numpy .zeros([n+1,2])
# u[:,0] = numpy.exp(-25*x**2)
# u[:,1] = u[:,0]
# u0 = numpy.zeros([n+1,2])
# uu = numpy.zeros([n+3,2])
# du = numpy.zeros([n+1, 2, nstage])
# for i in range(1001):
#   if i % 20 == 0:
#     print '%d / %d' % (i, 20000)
#     pylab.cla()
#     pylab.plot(x,u,'+-')
#     pylab.ylim([-2.0,2.0])
#     pylab.title('t = %f' % (i * dt))
#     pylab.savefig('wave_central_%04d.png' % int(i*dt*10/2))
#   u0[:,:] = u
#   for istage in range(nstage):
#     uu[1:-1,:] = u[:,:]
#     uu[0,:] = 2*u[0,:] - u[1,:]
#     uu[-1,:] = 2*u[-1,:] - u[-2,:]
#     du[:,0,istage] = (uu[2:,1] - uu[:-2,1]) / (2.0*dx)
#     du[:,1,istage] = (uu[2:,0] - uu[:-2,0]) / (2.0*dx)
#     u[:,:] = u0
#     for jstage in range(istage + 1):
#       u += du[:,:,jstage] * dt * stages.beta[istage][jstage]
#     u[0,0] = 0.0; u[-1,0] = 0.0
#     # u[0,1] = u[1,1]; u[-1,1] = u[-2,1]
# os.system('mencoder mf://wave_central_*.png -mf type=png:fps=10 -ovc lavc ' + \
#           '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o wave_central.avi')
# os.system('rm -rf wave_central_*.png')


# interpolation based derivative
gamma = 2.0 / dx
da1 = numpy.zeros([n+1, n+1])
da2 = numpy.zeros([n+1, n+1])
for i in range(n+1):
  if i < n - i:
    da1[:2*i+1,i] = interp_1d_grad_coef(i,x[:2*i+1],gamma)
    da2[:,i] = interp_1d_grad_coef(i,x,gamma)
  else:
    da1[:,i] = interp_1d_grad_coef(i,x,gamma)
    da2[2*i-n:,i] = interp_1d_grad_coef(n-i,x[2*i-n:],gamma)
#pylab.plot(da[:,0], '+-')
#pylab.plot(da[:,n/2], '+-')
#stop

u = numpy.zeros([n+1, 2])
u[:,0] = numpy.exp(-25*x**2)
u0 = numpy.zeros([n+1, 2])
du = numpy.zeros([n+1, 2, nstage])
for i in range(1001):
  if i % 10 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.ylim([-2.0,2.0])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('wave_interp_%04d.png' % int(i*dt*10))
  u0[:,:] = u
  for istage in range(nstage):
    for j in range(n+1):
      du[j,0,istage] = numpy.dot(da2[:,j], u[:,0])
      du[j,1,istage] = - numpy.dot(da1[:,j], u[:,1])
    u[:,:] = u0
    for jstage in range(istage + 1):
      u += du[:,:,jstage] * dt * stages.beta[istage][jstage]
    u[-1,0] = -u[-1,1]; u[0,1] = -u[0,0]
# os.system('mencoder mf://wave_interp_*.png -mf type=png:fps=10 -ovc lavc ' + \
#           '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o wave_interp.avi')
# os.system('rm -rf wave_interp_*.png')

