import pickle

import pylab
import numpy
from numpy import zeros, ones

from laplace import poisson

source, bcmask = pickle.loads(file('save.txt').read())
bcvalue = zeros(source.shape, dtype=float)

phi = poisson(source, bcmask, bcvalue)

pylab.contour(phi)
pylab.show()
