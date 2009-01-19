from math import fabs, log

import pylab

from neu import Neural


class MyNeural(Neural):
    target = 0.5
    d_target = 0.0011
    
    hist = ([], [])
    perf = []

    def monitor(self):
        self.target += self.d_target
        if self.target > 1.0:
            self.d_target = -fabs(self.d_target)
        elif self.target < 0.2:
            self.d_target = fabs(self.d_target)

        self.inputs[0].value = min(max(self.target, 0.0), 1.0)

        self.perf.append(self.performance())

    def performance(self):
        s = 0.0
        for n in self.nodes[:10]:
            s += n.output()
        self.hist[0].append(self.target)
        self.hist[1].append(s / 10)
        s = (1.1 - self.target) - s / 10
        if s**2 < 1E-8:
            return 0.0
        else:
            return -log(s**2)

    def hist_weights(self):
        weights = []
        for n in self.nodes:
            weights += n._weight
        pylab.hist(weights)

    def hist_output(self):
        outputs = []
        for n in self.nodes:
            outputs.append(n.output())
        pylab.hist(outputs, bins=100)

    def plot_output(self):
        i1, output1 = [], []
        i2, output2 = [], []
        for i, n in enumerate(self.nodes):
            if self.inputs[0] in n._input:
                i1.append(i)
                output1.append(n.output())
            else:
                i2.append(i)
                output2.append(n.output())
        pylab.plot(i1, output1, '+r')
        pylab.plot(i2, output2, '+b')
            



n = MyNeural(num_nodes=[30, 30, 30], num_connections=5, num_inputs = 1,
             num_input_duplication = 10)
for i in range(100):
    n.simulate(5000)
    # print n.target
    # pylab.cla()
    # n.plot_output()
    pylab.close('all')
    n.hist_weights()
    pylab.figure()
    n.plot_output()
    pylab.figure()
    pylab.plot(n.hist[0])
    pylab.plot(n.hist[1])
    pylab.ylim([min(min(n.hist[0]), min(n.hist[1])),
                max(max(n.hist[0]), max(n.hist[1]))])
    try: input()
    except: pass
