import numpy
import pylab

from numpy import asarray

from lorenz import Lorenz

l = Lorenz(10.0, 28.0, 8.0/3.0)
DT, NSTEP, NITER, STEP = 0.01, 1000, 100, 1.0E-11
ctrl = [0] * (NSTEP + 1)

obj0 = 1.0E+18
for iiter in range(NITER):
   traj = []
   x, y, z = 0.0, 10.0, 10.0
   
   traj.append((x, y, z))
   for i in range(NSTEP):
      x += ctrl[i]
      x, y, z = l.step(x, y, z, DT)
      traj.append((x, y, z))
   # obj = (asarray(traj)[int(NSTEP/2):,:]**2).sum()
   obj = (x+1)**2 + (y+1)**2 + (z-20)**2

   # if iiter % 100 == 0 or iiter == NITER - 1:
   #    t = asarray(traj[int(NSTEP/2):])
   #    pylab.figure()
   #    pylab.subplot(2,2,1)
   #    pylab.plot(t[:,0], t[:,1])
   #    pylab.xlabel('x'); pylab.ylabel('y')
   #    pylab.subplot(2,2,2)
   #    pylab.plot(t[:,0], t[:,2])
   #    pylab.xlabel('x'); pylab.ylabel('z')
   #    pylab.subplot(2,2,3)
   #    pylab.plot(t[:,1], t[:,2])
   #    pylab.xlabel('y'); pylab.ylabel('z')
   #    pylab.title('iteration %d' % iiter)
   
   if obj >= obj0:
      STEP = STEP / 2.0
      ctrl = [ctrl0[i] - trajadj[i][0] * STEP for i in range(NSTEP + 1)]
      continue
   else:
      STEP = STEP * 2.0
      obj0 = obj
      ctrl0 = ctrl
   print iiter, obj, STEP

   trajadj = []
   x, y, z = traj[-1]
   xadj, yadj, zadj = 2.0 * (x+1), 2.0 * (y+1), 2.0 * (z-20)
   trajadj.append((xadj, yadj, zadj))
   for i in range(NSTEP-1, 0-1, -1):
      x, y, z = traj[i]
      xadj, yadj, zadj = l.adj(x, y, z, xadj, yadj, zadj, DT)
      # if i >= int(NSTEP/2):
      #    xadj, yadj, zadj = xadj + 2.0 * x, yadj + 2.0 * y, zadj + 2.0 * z
      trajadj.append((xadj, yadj, zadj))
   trajadj.reverse()

   trajadj = numpy.asarray(trajadj)
   # trajadj[:,0] *= numpy.exp(numpy.arange(1001)/80)
   # trajadj[:,1] *= numpy.exp(numpy.arange(1001)/80)
   # trajadj[:,2] *= numpy.exp(numpy.arange(1001)/80)
   ctrl = [ctrl0[i] - trajadj[i][0] * STEP for i in range(NSTEP + 1)]
