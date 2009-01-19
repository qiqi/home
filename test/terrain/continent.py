from numpy import random, zeros, arange, exp, power, sqrt, absolute

def randimage(pertsize):
   a0 = random.random([3,3])
   for iter in range(len(pertsize)):
      a0[0,:] = 0.0
      a0[-1,:] = 0.0
      a0[:,0] = 0.0
      a0[:,-1] = 0.0
   
      newsize = [a0.shape[0]*2-1, a0.shape[1]*2-1]
   
      p = random.random(newsize) - 0.5
      # for i = 1:0
      #    p[2:2:end,:) = (p(1:2:end-1,:) + 2*p(2:2:end,:) + p(3:2:end,:)) / 4
   
      a = zeros(newsize)
      a[0::2, 0::2] = a0
      a[0::2, 1::2] = (a[0::2, 0:-1:2] + a[0::2, 2::2]) / 2.0
      a[1::2, 0::2] = (a[0:-1:2, 0::2] + a[2::2, 0::2]) / 2.0
      a[1::2, 1::2] = (a[1::2, 0:-1:2] + a[1::2, 2::2]) / 2.0
      a = a + p * 0.5**iter * pertsize[iter]
      a[0::2, 0::2] = a0
      a0 = a
   
   a = (a - a.min()) / (a.max() - a.min())
   return a

def continent(seed, itrans, sea, niter):
   # =========================================================
   #              Generate random continent
   # =========================================================
   # random.seed(seed=4); itrans, sea = 5, 0.465    # penesulas
   # random.seed(seed=7); itrans, sea = 7, 0.     # mountainous
   # random.seed(seed=13); itrans, sea = 13, 0.     # good plains
   # =========================================================

   pertsize = [1.0, 2.0, 4.0, 6.0, 4.0, 2.0, 1.0, 0.5, 0.2, 0.2, 0.1, 0.05]
   assert niter <= len(pertsize)

   random.seed(seed=seed)
   h = randimage(pertsize[0:niter])
   
   # =========================================================
   #              Generate transformation function
   # =========================================================
   from transform import transform
   h = transform(h,itrans)
   
   # =========================================================
   #              Flat out the sea
   # =========================================================
   if sea < h[0,0]:
      # find suitable sea level
      sea = h[1,1] + 0.01
   h = (h - sea) / (1 - sea)
   is_sea = (h < 0)
   h[is_sea] = 0.0

   # =========================================================
   #              Plateau variation
   # =========================================================
   pertsize = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 64.0, 64.0, 32.0, 16.0]
   gamma = randimage(pertsize[0:niter])
   b = exp((gamma - 0.3))
   h = power(h * 0.6 + 0.2, b) - power(0.2, b)
   h = h / h.max()
   h[is_sea] = 0.0
   # import pylab; pylab.contour(b); pylab.show()

   # =========================================================
   #              Add grid perturbation
   # =========================================================
   dh1, dh2 = zeros(h.shape), zeros(h.shape)
   dh1[:,0:-1] = absolute(h[:,0:-1] - h[:,1:])
   dh1[:,-1] = dh1[:,-2]
   dh1[:,1:-1] = (dh1[:,1:-1] + dh1[:,0:-2]) / 2
   dh2[0:-1,:] = absolute(h[0:-1,:] - h[1:,:])
   dh2[-1,:] = dh2[-2,:]
   dh2[1:-1,:] = (dh2[1:-1,:] + dh1[0:-2,:]) / 2
   dh = dh1 + dh2
   dh = (dh - dh.min()) / (dh.max() - dh.min())

   p = random.random(h.shape) * sqrt((dh * 0.9 + h * 0.1) * h)
   p[h == 0] = 0
   h += p * 2.0

   h /= h.max()
   h *= (1.0 - h.mean() * 0.7)
   
   return h

def slope(h):
   # calculate slope of terrain
   dh1, dh2 = zeros(h.shape), zeros(h.shape)
   dh1[:,0:-1] = absolute(h[:,0:-1] - h[:,1:])
   dh1[:,-1] = dh1[:,-2]
   dh1[:,1:-1] = (dh1[:,1:-1] + dh1[:,0:-2]) / 2
   dh2[0:-1,:] = absolute(h[0:-1,:] - h[1:,:])
   dh2[-1,:] = dh2[-2,:]
   dh2[1:-1,:] = (dh2[1:-1,:] + dh1[0:-2,:]) / 2
   dh = dh1 + dh2
   dh = (dh - dh.min()) / (dh.max() - dh.min())
   return dh
