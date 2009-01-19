import pylab

from continent import continent
from weather import temperature, precipitation

print 'loading...'

h = continent(13,13,0.0,8)    # good plains

# s, w = temperature(h*8, 10, 60)
# pylab.contour(h,[0])
# pylab.contour(s,[0],colors='r')
# pylab.axis('scaled')
# pylab.xlim([0, h.shape[0]-1])
# pylab.ylim([0, h.shape[0]-1])

p = precipitation(h, 10, 60)

print 'finished.'

pylab.contour(h,[0,0.001])
pylab.contour(p,[0.4,0.2,0.1,0.05,0.025,0.01])
pylab.axis('scaled')
pylab.xlim([0, h.shape[0]-1])
pylab.ylim([0, h.shape[0]-1])
