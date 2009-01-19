import os
import sys
import logging

import pylab

from leg import MonoLeg

m = [0.5,0.2,0.2]     # mass of the three pieces (kg)
l = [0.4,0.4,0.4]     # length of the three pieces (m)
i = [0.06,0.03,0.03]  # moment of inertia about the center (kg m^2)

pylab.matplotlib.interactive(False)
logging.root.setLevel(logging.DEBUG)

files =  []
leg = MonoLeg(m, l, i)
N = 100
DT = 0.1
for i in range(N):
   # sys.stdout.write('\x0d%d/%d' % (i, N))
   # sys.stdout.flush()
   if i * DT < 0.5:
      leg._control[0:2] = 0.1, -0.1
   else:
      leg._control[0:2] = 0.0

   leg.step(DT)
   print i, leg._a_pos[1] - leg._a_pos[2], leg._V[5] - leg._V[8]

   leg.plot([-3,3], [-2,4])
   fname = '_tmp%03d.png' % i
   pylab.savefig(fname)
   files.append(fname)

print
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc \
           -msglevel all=1 -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")
for fname in files:
   os.system("rm -f " + fname)

