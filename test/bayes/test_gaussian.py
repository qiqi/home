import pylab
import numpy
from numpy import random, dot, linalg, sqrt, exp, pi

from node import Node, Input

MEAN = numpy.array([2.0, 1.0])
COV = numpy.array([[2.0, 0.5], [0.5, 1.0]])

ndim = len(MEAN)
inputs = [Input() for i in range(ndim)]
output = Node(inputs)

errors = []
for i in range(1000):
    # sample a distribution
    vec = random.multivariate_normal(MEAN, COV)
    for x, input in zip(vec, inputs):
        input.value = x
    f = output.evaluate(i)
    # calculate error
    diff = vec - MEAN
    exponent = -0.5 * dot(diff, linalg.solve(COV, diff))
    prefix = 1.0 / sqrt((2*pi) ** ndim * linalg.det(COV))
    f0 = prefix * exp(exponent)
    errors.append(f - f0)

pylab.plot(errors)
