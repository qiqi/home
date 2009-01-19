import pylab
import numpy

randlist = []
for sty in range(10):
  fulllist = [numpy.random.random(200), numpy.random.random(200)]
  randlist.append([])
  for n in range(201):
    randlist[sty].append([fulllist[0][:n], fulllist[1][:n]])

def display(density, sty, animation):
  x, y = [], []
  for i in range(density.shape[0]):
    for j in range(density.shape[1]):
      n = int(round(density[i][j]))
      s = sty[i][j]
      x.extend((randlist[s][n][0] + animation[0]) % 1 + i)
      y.extend((randlist[s][n][1] + animation[1]) % 1 + j)
  pylab.cla()
  pylab.plot(x, y, '.', markersize=1)
  pylab.axis([0, density.shape[0], 0, density.shape[1]])

N = 50
pylab.matplotlib.interactive(False)
sty = numpy.random.randint(0,10,[N,N])
density = numpy.zeros([N,N], dtype=float);
spd = numpy.zeros([N,N], dtype=float);
density[5:15, 2:5] = 100.0;

speed = 0.1
for i in range(100):
  print i
  front = (density[:,1:] < density[:,:-1])
  spd = speed - front * (100 - density[:,1:]) * 0.0008
  density[:,1:] += density[:,:-1] * spd
  density[:,:-1] -= density[:,:-1] * spd
  display(density, sty, [0,i*speed])
  pylab.savefig('ani%02d.png' % i)

