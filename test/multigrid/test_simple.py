import pylab
import numpy
from numpy import zeros, ones

from laplace import poisson

# n = 128
n = 1024
# n = 4096
source = ones([n, n], dtype=float)
bcmask = zeros([n, n], dtype=bool)
bcvalue = zeros([n, n], dtype=float)
# bcmask[0,:] = True
# bcmask[-1,:] = True
# bcmask[:,0] = True
# bcmask[:,-1] = True

phi = poisson(source, bcmask, bcvalue)

# pylab.contour(phi)
# pylab.show()
