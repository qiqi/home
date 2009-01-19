import pylab
import numpy
from numpy import random, dot, linalg, sqrt, exp, cos, sin, pi

from node import Node, Input

MEAN = numpy.array([2.0, 1.0])
COV = numpy.array([[2.0, 0.5], [0.5, 1.0]])

def point_gen():
    a = random.uniform(0, 0.5 * pi)
    return numpy.array([cos(a), sin(a)])

# construct network
inputs = [Input(), Input()]
inter0 = Node([inputs[0]])
inter1 = Node([inputs[1]])
output = Node([inter0, inter1])

# train
for i in range(1000):
    # sample a distribution
    vec = point_gen()
    for x, input in zip(vec, inputs):
        input.value = x
    output.evaluate(i)

# evaluate performance
xvec = numpy.arange(-1, 2, 0.05)
yvec = numpy.arange(-1, 2, 0.05)
f = numpy.zeros((len(xvec), len(yvec)))
for ix, x in enumerate(xvec):
    for iy, y in enumerate(yvec):
        i += 1
        vec = numpy.array([x, y])
        for v, input in zip(vec, inputs):
            input.value = v
        f[ix, iy] = output.evaluate_test(i)

pylab.contour(xvec, yvec, f)
a = numpy.arange(0, 0.5 * pi, 0.01)
pylab.plot(cos(a), sin(a), '+b')
pylab.axis('equal')
