import sys

import pylab
import numpy
from numpy import zeros, arange, round

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from continent import continent, slope
from weather import temperature, precipitation
from vegetation import vegetation

# heightmap
h = continent(4,5,0.465,8)    # penesulas
# h = continent(7,7,0.0,8)      # mountainous
# h = continent(13,13,0.0,8)    # good plains
# h = continent(24,65,0.0,8)    # sahara
h *= 10.0

summer, winter = temperature(h, 10, 60)
precip = precipitation(h, 10, 60)
dh = slope(h)
r, g, b, t = vegetation(h, dh, summer, winter, precip)

color = zeros(h.shape + (3,))
color[:,:,0] = r
color[:,:,1] = g
color[:,:,2] = b

pylab.figure()
pylab.imshow(color)
pylab.axis('scaled')
pylab.xlim([0, h.shape[0]-1])
pylab.ylim([0, h.shape[0]-1])

pylab.figure()
h[h == 0] = -6.0
pylab.contour(h, list(arange(-6,0,1)) + \
                 list(arange(0.0,2.5,0.5)) + list(arange(2.5,10.0,1.0)))
pylab.axis('scaled')
pylab.xlim([0, h.shape[0]-1])
pylab.ylim([0, h.shape[0]-1])
# pylab.show()

h[h < 0] = 0.0
f = open('terrain.out', 'w')
for i in range(h.shape[0]):
   print i
   for j in range(h.shape[1]):
      f.write("%d %d %f %d %f %f %f\n" % ((i, j, h[i,j], t[i,j]) + \
                                          tuple(color[i,j,:])))
f.close()
