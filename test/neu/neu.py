from math import log
from random import sample

import pylab

from node import OrNode, AndNode, NotAndNode, Input



class Neural:
    nodes = []
    inputs = []

    def __init__(self, num_nodes, num_connections, num_inputs,
                 num_input_duplication = 1):
        num_or_nodes, num_and_nodes, num_notand_nodes = num_nodes
        for i in range(num_or_nodes):
            self.nodes.append(OrNode())
        for i in range(num_and_nodes):
            self.nodes.append(AndNode())
        for i in range(num_notand_nodes):
            self.nodes.append(NotAndNode())
        for i in range(num_inputs):
            self.inputs.append(Input())

        self.num_nodes = num_or_nodes + num_and_nodes + num_notand_nodes
        self.num_connections = num_connections
        self.num_inputs = num_inputs

        dup = [0] * self.num_inputs
        while min(dup) < num_input_duplication:
            print 'Generating connections with %d dups...' % \
                  num_input_duplication
            dup = [0] * self.num_inputs
            for i, n in enumerate(self.nodes):
                candidates = range(i) + range(i + 1, self.num_nodes) + \
                             range(self.num_nodes, self.num_nodes + \
                                   self.num_inputs * num_input_duplication)
                select = sample(candidates, num_connections)
                input = []
                for j in select:
                    if j < self.num_nodes:
                        input.append(self.nodes[j])
                    else:
                        j_input = (j - self.num_nodes) % self.num_inputs
                        input.append(self.inputs[j_input])
                        dup[j_input] += 1
                n.set_input(input)
        print 'Connections successfully generated with %d dups.' % min(dup)


    def compute(self):
        for node in self.nodes: node.compute()
        for node in self.nodes: node.update()

    def reward(self, amount):
        for node in self.nodes:
            node.reward(amount)

    def performance(self):
        return 0.0

    def monitor(self):
        pass

    def simulate(self, niter, inner=10):
        for i in xrange(niter):
            for j in xrange(inner):
                self.compute()
            self.reward(self.performance())
            self.monitor()
        

