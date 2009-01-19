import numpy
import pylab

from numpy import asarray

from lorenz import Lorenz

l = Lorenz(10.0, 28.0, 8.0/3.0)
DT, NSTEP = 0.01, 100

delta = 0.1
deltas = []
grads = []
for i in range(5):
   dtraj = []
   traj0 = []
   traj1 = []
   dx, dy, dz = 0.0, 10.0, 10.0
   x0, y0, z0 = 0.0, 10.0, 10.0
   x1, y1, z1 = x0 + delta * dx, y0 + delta * dy, z0 + delta * dz
   
   dtraj.append((dx, dy, dz))
   traj0.append((x0, y0, z0))
   traj1.append((x1, y1, z1))
   for i in range(NSTEP):
      dx, dy, dz = l.diff(x0, y0, z0, dx, dy, dz, DT)
      x0, y0, z0 = l.step(x0, y0, z0, DT)
      x1, y1, z1 = l.step(x1, y1, z1, DT)
      dtraj.append((dx, dy, dz))
      traj0.append((x0, y0, z0))
      traj1.append((x1, y1, z1))
   
   grad = (asarray(traj1[-1]) - asarray(traj0[-1])) / asarray(dtraj[-1])
   deltas.append(delta); grads.append(grad)
   delta *= 0.1

deltas, grads = asarray(deltas), asarray(grads)
pylab.loglog(deltas, grads[:,0],'+-')
pylab.loglog(deltas, grads[:,1],'+-')
pylab.loglog(deltas, grads[:,2],'+-')
