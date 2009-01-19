import copy

from math import log, exp

import random



ACTIVITY_DECAY = 0.9
ACTIVITY_EFFECT = 1.0E-5
ACTIVITY_RAND = 0.8E-1
OUTPUT_MIN = 1.0E-8
WEIGHT_MIN = 0.1

class Node:
    _output = 1.0
    _old_output = 1.0
    _input = []
    _weight = []
    _activity = []

    def reset(self):
        self._input = []

    def set_input(self, input_list):
        assert len(self._input) == 0
        self._input = copy.copy(input_list)
        self._weight = [1.0 / len(input_list)] * len(input_list)
        self._activity = [0.0] * len(input_list)

    def update(self):
        self._old_output = self._output
        for i in range(len(self._activity)):
            self._activity[i] *= ACTIVITY_DECAY

    def output(self):
        return self._output

    def reward(self, amount):
        total_weight = 0.0
        for i in range(len(self._weight)):
            try:
                self._weight[i] *= exp(random.gauss(0, ACTIVITY_RAND) + \
                                   amount * self._activity[i] * ACTIVITY_EFFECT)
            except OverflowError:
                print '!--'
                print random.gauss(0, ACTIVITY_RAND), amount, self._activity[i]
                print '--!'
            self._weight[i] = max(self._weight[i], WEIGHT_MIN)
            total_weight += self._weight[i]
        # normalize weights
        for i in range(len(self._weight)):
            self._weight[i] /= total_weight
        

    def compute(self):
        assert len(self._input) > 0
        result = self.do_compute()
        result = max(result, OUTPUT_MIN)
        self._output = result
        


class OrNode(Node):
    def do_compute(self):
        result = 0.0
        for i in range(len(self._input)):
            a = self._input[i].output() * self._weight[i]
            if self._input[i].output() > 1 or self._input[i].output() < 0:
                print self._input[i], self._input[i].output(), \
                      [j.output() for j in self._input[i]._input]
            self._activity[i] += a
            result += a
        return result



class AndNode(Node):
    def do_compute(self):
        result = 0.0
        for i in range(len(self._input)):
            a = log(self._input[i].output()) * self._weight[i]
            self._activity[i] += a**2
            result += a
        return result



class NotAndNode(Node):
    def do_compute(self):
        result = 0.0
        for i in range(len(self._input)):
            a = log(self._input[i].output()) * self._weight[i]
            self._activity[i] += a**2
            result = a
        result = exp(result)
        return 1.0 - result



class Input:
    value = 0.5
    def output(self):
        return self.value

