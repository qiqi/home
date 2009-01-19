import numpy
import pylab
from numpy import zeros, ones, linspace, kron, linalg, exp, sqrt, diag, \
                  arctan, pi

from arbinterp import interp_nd_coef, max_order_set, boundary_set_and_zeta

pylab.matplotlib.interactive(False)

gamma = [1, 10, 25, 1000]
nz = 25; z = linspace(0, 1, nz)

xall = []
for line in file('quasi-random/niederreiter_2d.txt').readlines():
  xall.append(map(float, line.strip().split()))
xall = numpy.array(xall)

cc = linspace(-2,3,51)
lw = ones(51); lw[[20,30]] = 5; lw = tuple(lw)

for n in numpy.arange(3, 4, 2)**2:
  print 'n = ', n
  x = xall[:n,:]
  for i in range(len(gamma)):
    print '    gamma = ', gamma[i]
    # calculate order set
    k = 0; order_set = []
    while len(order_set) < n:
      k += 1
      order_set = max_order_set(2, k)
    boundary_set, boundary_zeta = boundary_set_and_zeta(order_set)
    # calculate basis
    a = zeros([nz,nz,n]); er2 = zeros([nz,nz])
    for jx in range(nz):
      for jy in range(nz):
        a[jx,jy,:], b, er2[jx,jy] = \
        interp_nd_coef(numpy.array([z[jx], z[jy]]), x, zeros(n), \
                       zeros(0), zeros(0), 1.0, gamma[i], \
                       order_set, boundary_set, boundary_zeta)
    # plot basis
    for j in range(n):
      pylab.figure(figsize=(10,10))
      pylab.contour(z, z, a[:,:,j], cc, linewidths=lw)
      for k in range(n):
        pylab.plot(x[:,1], x[:,0], 'ok', mfc='w')
      pylab.plot([x[j,1]], [x[j,0]], 'ok')
      pylab.axis('scaled'); pylab.axis([-0.1,1.1,-0.1,1.1])
      pylab.savefig('basis_2d/random_%07.1fn%03d_%03d.png' % (gamma[i], n, j))


