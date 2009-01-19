import numpy
import pylab

from numpy import asarray

from lorenz import Lorenz

l = Lorenz(10.0, 28.0, 8.0/3.0)
dtraj, traj = [], []
DT, NSTEP = 0.01, 1000
dx, dy, dz = 0.0, 10.0, 10.0
x, y, z = 0.0, 10.0, 10.0

dtraj.append((dx, dy, dz))
traj.append((x, y, z))
for i in range(NSTEP):
   dx, dy, dz = l.diff(x, y, z, dx, dy, dz, DT)
   x, y, z = l.step(x, y, z, DT)
   dtraj.append((dx, dy, dz))
   traj.append((x, y, z))

trajadj = []
xadj, yadj, zadj = 0.0, 10.0, 10.0
trajadj.append((xadj, yadj, zadj))
for i in range(NSTEP-1, 0-1, -1):
   x, y, z = traj[i]
   xadj, yadj, zadj = l.adj(x, y, z, xadj, yadj, zadj, DT)
   trajadj.append((xadj, yadj, zadj))
trajadj.reverse()

aggr = (asarray(trajadj) * asarray(dtraj)).sum(axis=1)
# pylab.plot(aggr)
print aggr.max() - aggr.min()

pylab.plot(numpy.arange(len(trajadj)) * DT, trajadj)
