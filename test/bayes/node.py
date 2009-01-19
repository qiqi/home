import numpy
from numpy import zeros, eye, asmatrix, dot, pi, exp, sqrt, linalg



class Input(object):

    def __init__(self):
        self.value = 0.0

    def evaluate(self, ieval):
        return self.value



class Node(object):

    def __init__(self, inputs):
        self._inputs = inputs
        self._ndim = len(inputs)
        self._ieval = -1
        self._result = 0
        # bayesian priors
        self._mu = zeros(self._ndim)
        self._sigma = asmatrix(eye(self._ndim))
        self._n = 0.1

    def evaluate(self, ieval):
        assert ieval >= self._ieval
        if self._ieval < ieval:
            self._ieval = ieval
            inputs = [node.evaluate(ieval) for node in self._inputs]
            inputs = numpy.array(inputs)
            # calculate posterior
            diff = inputs - self._mu
            self._mu = (self._n * self._mu + inputs) / (self._n + 1.0)
            sigma_new = self._n / (self._n + 1.0) * \
                        asmatrix(diff).transpose() * asmatrix(diff)
            self._sigma = (self._n * self._sigma + sigma_new) / (self._n + 1.0)
            self._n += 1.0 - self._n * 1.0E-5
            # calculate result
            diff = inputs - self._mu
            exponent = -0.5 * dot(diff, linalg.solve(self._sigma, diff))
            prefix = 1.0 / sqrt((2*pi) ** self._ndim * linalg.det(self._sigma))
            self._result = prefix * exp(exponent)
        return self._result

    def evaluate_test(self, ieval):
        assert ieval >= self._ieval
        if self._ieval < ieval:
            self._ieval = ieval
            inputs = [node.evaluate(ieval) for node in self._inputs]
            inputs = numpy.array(inputs)
            # calculate result
            diff = inputs - self._mu
            exponent = -0.5 * dot(diff, linalg.solve(self._sigma, diff))
            prefix = 1.0 / sqrt((2*pi) ** self._ndim * linalg.det(self._sigma))
            self._result = prefix * exp(exponent)
        return self._result

