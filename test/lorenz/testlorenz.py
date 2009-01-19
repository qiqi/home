import numpy
import pylab

from lorenz import Lorenz

l = Lorenz(10.0, 28.0, 8.0/3.0)
traj = []
DT, NSTEP = 0.01, 10000
x, y, z = 0.0, 10.0, 10.0

traj.append((x, y, z))
for i in range(NSTEP):
   x, y, z = l.step(x, y, z, DT)
   traj.append((x, y, z))

traj = numpy.asarray(traj)
pylab.subplot(2,2,1)
pylab.plot(traj[:,0], traj[:,1])
pylab.xlabel('x'); pylab.ylabel('y')
pylab.subplot(2,2,2)
pylab.plot(traj[:,0], traj[:,2])
pylab.xlabel('x'); pylab.ylabel('z')
pylab.subplot(2,2,3)
pylab.plot(traj[:,1], traj[:,2])
pylab.xlabel('y'); pylab.ylabel('z')

