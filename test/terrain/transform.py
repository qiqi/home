from numpy import random, arange, exp, zeros

def transform(a, icase):
   # Generate transformation function
   random.seed(icase)

   njump = 4        # big jumps
   pos = random.random(size=njump)
   mag = random.random(size=njump)
   shp = random.random(size=njump)

   a = (a - a.min()) / (a.max() - a.min())
   b = zeros(a.shape)
   x = arange(0, 1.001, 0.001)
   y = zeros(x.shape)
   for i in range(njump):
      e = exp((a - pos[i]) * 50 * shp[i])
      b = b + mag[i] * (e - 1/e) / (e + 1/e)
      e = exp((x - pos[i]) * 50 * shp[i])
      y = y + mag[i] * (e - 1/e) / (e + 1/e)
   b = (b - y[0]) / (y[-1] - y[0])
   y = (y - y[0]) / (y[-1] - y[0])

   nsmall = 10    # small jumps
   pos = random.random(size=nsmall)
   mag = random.random(size=nsmall) / 10.0
   shp = random.random(size=nsmall) * 25.0
   for i in range(nsmall):
      imin = (pos[i] < y).sum()
      pos[i] = x[imin]

   for i in range(nsmall):
      e = exp((a - pos[i]) * 25 * shp[i])
      b = b + mag[i] * (e - 1/e) / (e + 1/e)
      e = exp((x - pos[i]) * 25 * shp[i])
      y = y + mag[i] * (e - 1/e) / (e + 1/e)
   a = (b - y[0]) / (y[-1] - y[0])
   y = (y - y[0]) / (y[-1] - y[0])

   #import pylab; pylab.figure(); pylab.plot(x,y)

   return (a - a.min()) / (a.max() - a.min())

