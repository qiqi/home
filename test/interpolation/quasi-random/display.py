import pylab

pylab.matplotlib.interactive(False)

x, y = [], []
for line in file('niederreiter_2d.txt').readlines():
  tx, ty = line.strip().split()
  x.append(float(tx))
  y.append(float(ty))

for i in range(16,5000,16):
  pylab.plot(x[:i], y[:i], 'ok')
  pylab.axis([-0.1, 1.1, -0.1, 1.1])
  pylab.savefig('n2d%05d.png' % i)
