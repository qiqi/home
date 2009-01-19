import os

import numpy
import pylab
from numpy import linspace

import rungekutta
from arbinterp import interp_1d_grad_coef

def niederreiter_1d(x0, x1, n):
  niederreiter = []
  for line in file('niederreiter_1d.txt', 'r').readlines()[:n]:
    niederreiter.append(float(line.strip()))
  return numpy.array(sorted(niederreiter)) * (x1 - x0) + x0

n = 50
x = niederreiter_1d(-1, 1, n+1)
xx = numpy.zeros(3*n+1)
xx[:n+1] = x + x[0] - x[-1]
xx[n:2*n+1] = x
xx[2*n:] = x + x[-1] - x[0]

dt = 0.01
stages = rungekutta.runge6
nstage = len(stages.alpha)

pylab.figure(figsize=(5,5))

# second order central difference
# dxp = xx[n+1:2*n+2] - xx[n:2*n+1]
# dxm = xx[n:2*n+1] - xx[n-1:2*n]
# cp = dxm / dxp / (dxm + dxp)
# cm = - dxp / dxm / (dxm + dxp)
# cc = - (cm + cp)
# 
# u = numpy.exp(-25*x**2)
# u0 = numpy.zeros(n+1)
# uu = numpy.zeros(3*n+1)
# du = numpy.zeros([n+1, nstage])
# for i in range(100001):
#   if i % 10000 == 0:
#     print '%d / %d' % (i, 100000)
#     pylab.cla()
#     pylab.plot(x,u,'+-')
#     pylab.plot(linspace(-1,1,1000),numpy.exp(-25*linspace(-1,1,1000)**2),':k')
#     pylab.ylim([-0.2,1.1])
#     pylab.title('t = %f' % (i * dt))
#     pylab.savefig('centralrand_%04d.png' % int(i*dt*10/2))
#   u0[:] = u
#   for istage in range(nstage):
#     uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
#     du[:,istage] = -cp * uu[n+1:2*n+2] - cc * uu[n:2*n+1] - cm * uu[n-1:2*n]
#     u[:] = u0
#     for jstage in range(istage + 1):
#       u += du[:,jstage] * dt * stages.beta[istage][jstage]
# os.system('mencoder mf://centralrand_*.png -mf type=png:fps=10 -ovc lavc ' + \
#           '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o centralrand.avi')
# os.system('rm -rf centralrand_*.png')
# stop


# interpolation based derivative
gamma = 2.0 / (x[1:] - x[:-1]).mean()
da = numpy.zeros([2*n+1, n+1])
for i in range(n+1):
  da[:,i] = interp_1d_grad_coef(n,xx[i:2*n+1+i],gamma)

u = numpy.exp(-25*x**2)
#u[:] = 1.0
#u[x>0] = 0.0
#u[-1] = u[0]
u0 = numpy.zeros(n+1)
uu = numpy.zeros(3*n+1)
du = numpy.zeros([n+1, nstage])
for i in range(20001):
  if i % 2000 == 0:
    print '%d / %d' % (i, 20000)
    pylab.cla()
    pylab.plot(x,u,'+-')
    pylab.plot(linspace(-1,1,1000),numpy.exp(-25*linspace(-1,1,1000)**2),':k')
    pylab.ylim([-0.2,1.1])
    pylab.title('t = %f' % (i * dt))
    pylab.savefig('interprand_%04d.png' % int(i*dt*10/2))
  u0[:] = u
  for istage in range(nstage):
    uu[:n] = u[:-1]; uu[n:2*n] = u[:-1]; uu[2*n:] = u
    for j in range(n+1):
      du[j,istage] = -numpy.dot(da[:,j], uu[j:2*n+1+j])
    u[:] = u0
    for jstage in range(istage + 1):
      u += du[:,jstage] * dt * stages.beta[istage][jstage]
os.system('mencoder mf://interprand_*.png -mf type=png:fps=10 -ovc lavc ' + \
          '-lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o interprand.avi')
# os.system('rm -rf interprand_*.png')

