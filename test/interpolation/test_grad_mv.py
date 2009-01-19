from numpy import *
from interp1d import *
from interp2d import *

order = 5
order_set = max_order_set(2, order)
boundary_set, boundary_zeta = boundary_set_and_zeta(order_set)
gamma = 7.0

for n in range(3,19,2):
  x = zeros([n*n,2])
  for i in range(n*n):
    x[i,0] = float(i/n - n/2) / (n-1)
    x[i,1] = float(i%n - n/2) / (n-1)
  ix = n/2 * (n+1)
  da = interp_nd_grad_coef(ix, x, gamma, order_set, boundary_set, boundary_zeta)
  
  x0 = [1,0.5]
  u = 1. / ((x[:,0] - x0[0])**2 + (x[:,1] - x0[1])**2)
  xi = x[ix,:]
  ux = -2 * (xi[0] - x0[0]) / ((xi[0] - x0[0])**2 + (xi[1] - x0[1])**2)**2
  uy = -2 * (xi[1] - x0[1]) / ((xi[0] - x0[0])**2 + (xi[1] - x0[1])**2)**2
  print n, ux, uy, dot(u, da)
