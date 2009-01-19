import copy
import sys
import pickle
import time

import numpy
import pylab
from numpy import zeros, ones, eye, kron, linalg, dot, exp, sqrt, diag, pi, \
                  asarray, sign



# -------------------------------------------------------------- #
#                                                                #
#                             1D STUFF                           #
#                                                                #
# -------------------------------------------------------------- #

def factorial(n):
  """
  Factorial of an integer or an integer list.
  For a list of integer, factorial([a,b,c]) = a! b! c!
  """
  if isinstance(n, int) or isinstance(n, float):
    if n <= 1:
      return 1.0
    else:
      return n * factorial(n-1)
  else:
    res = 1.0
    for ni in n:
      res *= factorial(ni)
    return res



def solve_L(L, b):
  """
  Solve lower triangular system
  """
  n = b.size
  assert L.shape == (n,n)
  x = zeros(n)
  for i in range(n):
    x[i] = (b[i] - dot(x[:i], L[i,:i])) / L[i,i]
    if not numpy.isfinite(x[i]):
      x[i] = 0.0
  return x



def solve_R(R, b):
  """
  Solve upper triangular system
  """
  n = b.size
  assert R.shape == (n,n)
  x = zeros(n, dtype=R.dtype)
  for i in range(n-1,-1,-1):
    x[i] = (b[i] - dot(x[i+1:], R[i,i+1:])) / R[i,i]
    if not numpy.isfinite(x[i]):
      x[i] = 0.0
  return x



def quad_prog(X, G, H, c):
  """
  Solve the constraint quadratic programming with
  min ab^T A ab, where A = X^T X + G^2 + H^2, s.t. ab^T c = 1
  To get best performance, order columns of X, G and H such that the diagonal
  of A is of increasing size.
  """
  # construct A
  A = dot(X.transpose(), X) + diag(G**2 + H**2)
  # normalize and filter out infinite columns and rows
  D = sqrt(diag(A))
  finite = numpy.isfinite(D)
  X = X[:,finite]
  H = H[finite]
  G = G[finite]
  # QR decomposition of X, R is now cholesky factor of A
  n1, n2 = X.shape
  A = zeros([n1+n2, n2])
  A[:n1,:] = X
  A[range(n1,n1+n2), range(n2)] = sqrt(G**2 + H**2)
  R = numpy.linalg.qr(A, mode='r')
  # solve...
  ab = zeros(c.size, dtype=X.dtype)
  tmp = solve_L(R.transpose(), c[finite])
  ab[finite] = solve_R(R, tmp)
  ab /= dot(c, ab)
  if numpy.isnan(ab).any():
    stop
  return asarray(ab)



def interp_1d_coef(z, x, dx, y, dy, beta, gamma):
  """
  Calculate interpolation coefficients in 1D.
  z is the point where the interpolation is evaluated.
  x is the points where the function value is available.
  dx is the measurement error at x, must has same size as x.
  y is the points where the function gradient is available.
  dy is the measurement error at y, must has same size as y.
  beta is the `magnitude' of the target function.
  gamma is the `wave number' of the target function.
    combined with beta, it provides an estimate of the derivative growth:
      f^(k) = O(beta * gamma**k)
    larger gamma = more conservative and lower order interpolation.
  return values:
    a, b, er2 = interp_1d_coef(...)
  a and b are the interpolation coefficients.
  er2 is the expected squared residual.
  """
  # verifying arguments
  gamma = float(gamma); z = float(z)
  assert x.ndim == dx.ndim == 1 and x.size == dx.size
  assert y.ndim == dy.ndim == 1 and y.size == dy.size
  n = x.size; m = y.size; N = min(n+m, 200)
  # calculate interpolant at z
  if abs(x - z).min() < 1.0E-12: # exactly matches a data point
    imatch = abs(x - z).argmin()
    if dx[imatch] == 0.0:        # and function value on that point is exact
      a = numpy.array([0]*imatch+[1]+[0]*(n-imatch-1))
      b = numpy.zeros(m)
      return a, b, 0.0
  # construct X = [Xa,Xb]
  X = zeros([N, n+m], dtype=float)
  for i in range(N):
    X[i,:n] = gamma**(i+1) / factorial(i+1) * (x - z) ** (i+1)
    X[i,n:] = gamma**(i+1) / factorial(i) * (y - z) ** i
  X *= beta
  # construct diagonal G matrix for the Lagrange residual
  G = zeros(n+m)
  G[:n] = gamma**(N+1) / factorial(N+1) * (x - z)**(N+1)
  G[n:] = gamma**(N+1) / factorial(N) * (y - z)**N
  G *= beta
  # construct diagonal H matrix for measurement errors
  H = zeros(n+m)
  H[:n] = dx
  H[n:] = dy
  # construct c
  c = zeros(n+m); c[:n] = 1.0
  # first try to assemble the diagonal of matrix A, and sort by its diagonal
  diagA = (X**2).sum(0) + G**2 + H**2
  isort = sorted(range(n+m), key=diagA.__getitem__)
  irevt = sorted(range(n+m), key=isort.__getitem__)
  # permute columns of X and diagonal of G and H
  X = X[:,isort]
  G = G[isort]
  H = H[isort]
  c = c[isort]
  ab = quad_prog(X, G, H, c)
  # reverse sorting permutation to get a and b and normalize
  abrevt = ab[irevt]
  a = abrevt[:n]
  b = abrevt[n:]
  # compute the expeted squared residual
  finite = (ab != 0)
  Xab = dot(X[:,finite], ab[finite])**2
  Gab = (G*ab)[finite]**2
  Hab = (H*ab)[finite]**2
  er2 = Xab.sum() + Gab.sum() + Hab.sum()
  return a, b, er2



def interp_1d_grad_coef(ix, x, gamma):
  """
  Calculate interpolant gradient coefficients in 1D.
  z is the point where the interpolant gradient is evaluated.
  x is the points where the function value is available.
  gamma is the `wave number' of the target function.
    it provides an estimate of the derivative growth:
      f^(k) = O(gamma**k)
    larger gamma = more conservative and lower order interpolation.
  return values:
    da = interp_1d_grad_coef(...)
  da contains the interpolant gradient coefficients.
  """
  # verifying arguments
  gamma = float(gamma)
  assert x.ndim == 1
  n = x.size; N = min(n, 50)
  # construct X
  X = zeros([N, n], dtype=float)
  for i in range(N):
    X[i,:] = gamma**(i+1) / factorial(i+1) * (x - x[ix]) ** (i+1)
  # construct diagonal G matrix for the Lagrange residual
  G = gamma**(N+1) / factorial(N+1) * (x - x[ix])**(N+1)
  # construct c
  c = ones(n)
  # first try to assemble the diagonal of matrix A, and sort by its diagonal
  diagA = (X**2).sum(0) + G**2
  finite = numpy.isfinite(diagA)
  isort = sorted(range(n), key=diagA.__getitem__)
  irevt = sorted(range(n), key=isort.__getitem__)
  # drop first one -- all 0
  assert isort[0] == ix
  isort = isort[1:]
  # permute columns of X and diagonal of G
  X = X[:,isort]
  G = G[isort]
  c = c[isort]
  finite = finite[isort]
  # solve linear system
  X = X[:,finite]
  G = G[finite]
  # QR decomposition of X, R is now cholesky factor of A
  n1, n2 = X.shape
  A = zeros([n1+n2, n2])
  A[:n1,:] = X
  A[range(n1,n1+n2), range(n2)] = G
  R = numpy.linalg.qr(A, mode='r')
  # ------------------- gradient part -----------------
  # calculate rhs = dA * a
  rhs = - gamma * X[0,:];
  # solve for da
  da = zeros(c.size + 1)
  tmp = solve_L(R.transpose(), rhs)
  da[1:][finite] = solve_R(R, tmp)
  da[0] = -dot(da[1:], c)
  # reverse sorting permutation to get a and b and normalize
  da = da[irevt]
  return -da



def calc_beta(fx, dfx):
  """
  Estimate the `magnitude' parameter beta from data points.
  """
  assert fx.ndim == 1 and fx.shape == dfx.shape
  n = fx.size
  f_bar = fx.mean()
  ratio = (dfx**2).sum() / ((fx - f_bar)**2).sum() * (n-1) / float(n)
  beta = sqrt(((fx - f_bar)**2).sum() / (n-1) * exp(-ratio))
  return beta



def calc_res_ratio_avg_1d(beta, gamma, x, fx, dfx, y, fpy, dfpy):
  """
  A utility function used by calc_gamma, calculates the average ratio of real
  residual to the estimated residual, which is used to make decision in the
  bisection.
  """
  n = x.size
  n_remove = max(1, n / 5)
  i_sorted = sorted(range(n), key=x.__getitem__)
  var_total = 0.0
  for i in range(0,n,n_remove):
    # interpolate each value data point using all other data points
    target = i_sorted[i:i+n_remove]
    base = i_sorted[:i] + i_sorted[i+n_remove:]
    for j in target:
      a, b, er2 = interp_1d_coef(x[j], x[base], dfx[base], y, dfpy, beta, gamma)
      # average the ratio of real residual to estimated residual
      resid = dot(a,fx[base]) + dot(b,fpy) - fx[j]
      var_total += resid**2 / (er2 + dfx[j]**2)
  return var_total / n



def calc_gamma_1d(x, fx, dfx, y, fpy, dfpy):
  """
  Estimate the `wave number' parameter gamma from data points.
  This function prints stuff.
  """
  beta = calc_beta(fx, dfx)
  print '    using beta = ', beta
  # upper and lower bounds
  sorted_x = numpy.array(sorted(x))
  sorted_y = numpy.array(sorted(y))
  if y.size == 0:
    delta_max = x.max() - x.min()
    delta_min = (sorted_x[1:] - sorted_x[:-1]).min()
  else:
    delta_max = max(x.max(), y.max()) - min(x.min(), y.min())
    delta_min = min((sorted_x[1:] - sorted_x[:-1]).min(), \
                    (sorted_y[1:] - sorted_y[:-1]).min())
  assert delta_max > delta_min
  gamma_min = 1. / delta_max
  gamma_max = pi / delta_min
  # logorithmic bisection for gamma
  while gamma_max / gamma_min > 1.1:
    print '    bisecting [', gamma_min, ',', gamma_max, '] for gamma...'
    gamma_mid = sqrt(gamma_max * gamma_min)
    res_ratio = calc_res_ratio_avg_1d(beta, gamma_mid, x, fx, dfx, y, fpy, dfpy)
    if res_ratio < 1.0:
      gamma_max = gamma_mid
    else:
      gamma_min = gamma_mid
  # final selected gamma
  gamma_mid = sqrt(gamma_max * gamma_min)
  print '    using gamma = ', gamma_mid
  return gamma_mid



def calc_gamma_1d_fast(x, fx, dfx, y, fpy, dfpy):
  """
  Estimate the `wave number' parameter gamma from data points.
  This function prints stuff.
  """
  n, m = x.size, y.size
  # sort x
  isort = sorted(range(x.size), key=x.__getitem__)
  x, fx, dfx = x[isort], fx[isort], dfx[isort]
  # calculate derivative from fx
  fpx = (fx[1:] - fx[:-1]) / (x[1:] - x[:-1])
  dfpx = sqrt(dfx[1:]**2 + dfx[:-1]**2) / (x[1:] - x[:-1])
  # estimate beta*gamma
  ratio = ((dfpx**2).sum() + (dfpy**2).sum()) / \
          ((fpx**2).sum() + (fpy**2).sum())
  betagamma = (((fpx**2).sum() + (fpy**2).sum()) / (n+m) * exp(-ratio)) ** 0.5
  # detect discontinuity, raise gamma if needed
  if m > 0:
    dfmax = max(abs(fpx).max(), abs(fpy).max())
  else:
    dfmax = abs(fpx).max()
  betagamma = max(betagamma, 0.16 * dfmax)
  gamma = betagamma / calc_beta(fx, dfx) * 4
  print '    using gamma = ', gamma
  return gamma



def interp_1d(z, x, fx, dfx=None, y=None, fpy=None, dfpy=None, \
              compute_dfz=False):
  """
  Wang interpolation in 1D.
  z is the points where the interpolation is evaluated.
  x is the points where the function value is available and given by fx, dfx
    estimates the error in the function values, default to 0 (exact value).
  y is the points where the function gradient is available and given by fpy,
    dfpy estimates the error in the gradient values, default to 0
    (exact gradient).
  compute_dfz indicates whether an estimate of error in the interpolation
    approximation is returned.  If False (default), "fz = interp_1d(...)";
    if true, "fz, dfz = interp_1d(...)"
  """
  # verifying and handling arguments
  assert z.ndim == x.ndim == fx.ndim == 1 and x.size == fx.size
  if dfx is None:
    dfx = zeros(x.shape)
  else:
    assert dfx.shape == x.shape
  if y is None:
    assert fpy is None and dfpy is None
    y = zeros(0); fpy = zeros(0); dfpy = zeros(0)
  else:
    assert fpy is not None and fpy.shape == y.shape
    if dfpy is None:
      dfpy = zeros(y.shape)
    else:
      assert dfpy.shape == y.shape
  n = x.size; m = y.size
  # calculate beta and gamma
  beta = calc_beta(fx, dfx)
  gamma = calc_gamma_1d(x, fx, dfx, y, fpy, dfpy)
  # interpolation for each z[i]
  t0 = time.time()
  fz, dz = [], []
  for zi in z:
    a, b, er2 = interp_1d_coef(zi, x, dfx, y, dfpy, beta, gamma)
    fz.append(dot(a, fx) + dot(b, fpy))
    dz.append(sqrt(er2))
  print '    time: ', time.time() - t0
  if compute_dfz:
    return numpy.array(fz), numpy.array(dz)
  else:
    return numpy.array(fz)



# -------------------------------------------------------------- #
#                                                                #
#                             2D STUFF                           #
#                                                                #
# -------------------------------------------------------------- #

def max_order_set(d, max_order):
  """
  Return an d-dimensional expansion order set S = {k | |k|<max_order}.
  """
  order_set = [(0,) * d]
  ptr = 0
  while ptr < len(order_set):
    kappa = order_set[ptr]
    if sum(kappa) < max_order:
      for i in range(d):
        kappa_new = kappa[:i] + (kappa[i] + 1,) + kappa[i+1:]
        if kappa_new not in order_set:
          order_set.append(kappa_new)
    ptr += 1
  return order_set



def boundary_set_and_zeta(order_set):
  """
  Return the boundary of an expansion order set.
  """
  d = len(order_set[0])
  boundary_set = []
  for kappa in order_set:
    for k in range(d):
      kappa_new = kappa[:k] + (kappa[k] + 1,) + kappa[k+1:]
      if kappa_new not in order_set and kappa_new not in boundary_set:
        boundary_set.append(kappa_new)
  boundary_zeta = []
  for kappa in boundary_set:
    zeta = ()
    for k in range(d):
      if kappa[k] > 0 and kappa[:k] + (kappa[k]-1,) + kappa[k+1:] in order_set:
        zeta += (kappa[k],)
      else:
        zeta += (0,)
    boundary_zeta.append(zeta)
  return boundary_set, boundary_zeta



def interp_nd_coef(z, x, dx, y, dy, beta, gamma, \
                   order_set, boundary_set, boundary_zeta):
  """
  Calculate interpolation coefficients in n-dimensional spacee.
  z is the point where the interpolation is evaluated, the size of z
    is the dimension of the space.
  x is the points where the function value is available, the columns of x
    must be the same as the size of z.
  dx is the measurement error at x, its size must be same as the rows of x.
  y is the points where the function gradient is available, the columns of y
    must be the same as the size of z.
  dy is the measurement error at y, must has same shape as y.
  beta is the `magnitude' of the target function.
  gamma is the `wave number' of the target function.
    combined with beta, it provides an estimate of the derivative growth:
      f^(k) = O(beta * gamma**k)
    larger gamma = more conservative and lower order interpolation.
  order_set is the expansion order set of the Taylor expansions.
  boundary_set is the boundary of order_set.
  boundary_zeta.
  return values:
    a, b, er2 = interp_nd_coef(...)
  a and b are the interpolation coefficients.
  er2 is the expected squared residual.
  """
  # verifying arguments
  gamma = float(gamma)
  d, n, m = z.size, x.shape[0], y.shape[0]
  assert z.shape == (d,)
  assert x.shape == (n, d) and dx.shape == (n,)
  assert y.shape == (m, d) and dy.shape == (m, d) or m == 0
  assert len(order_set) > 1 and order_set[0] == (0,) * d
  # calculate interpolant at z
  d2 = ((x - z)**2).sum(1)
  if d2.min() < 1.0E-18:    # exactly matches a data point
    imatch = d2.argmin()
    if dx[imatch] == 0.0:   # and function value on that point is exact
      a = numpy.array([0]*imatch+[1]+[0]*(n-imatch-1))
      b = numpy.zeros([m, d])
      return a, b, 0.0
  # construct X = [Xa,Xb]
  N = len(order_set) - 1
  X = zeros([N, n+m*d], dtype=float)
  for i, kappa in enumerate(order_set[1:]):
    # Xa part
    X[i,:n] = gamma**sum(kappa) / factorial(kappa) * ((x - z)**kappa).prod(1)
    # Xb part
    if m > 0:
      for k in range(d):
        if kappa[k] > 0:
          kappa_p = kappa[:k] + (kappa[k] - 1,) + kappa[k+1:]
          X[i,n+k*m:n+(k+1)*m] = gamma**sum(kappa) / factorial(kappa_p) * \
                                 ((y - z)**kappa_p).prod(1)
  X *= beta
  # construct diagonal G matrix for the Lagrange residual
  G2 = zeros(n+m*d)
  for kappa, zeta in zip(boundary_set, boundary_zeta):
    Gi = gamma**sum(kappa) / factorial(kappa) * \
         ((x - z)**kappa).prod(1) * sum(zeta) / sum(kappa)
    G2[:n] += Gi**2
    if m > 0:
      for k in range(d):
        if kappa[k] > 0:
          kappa_p = kappa[:k] + (kappa[k] - 1,) + kappa[k+1:]
          zeta_p = zeta[:k] + (max(0, zeta[k] - 1),) + zeta[k+1:]
          Gi = gamma**sum(kappa) / factorial(kappa_p) * \
               ((y - z)**kappa_p).prod(1) * sum(zeta_p) / sum(kappa_p)
          G2[n+k*m:n+(k+1)*m] += Gi**2
  G = sqrt(G2) * beta
  # construct diagonal H matrix for measurement errors
  H = zeros(n+m*d)
  H[:n] = dx
  if m > 0:
    for k in range(d):
      H[n+k*m:n+(k+1)*m] = dy[:,k]
  # construct c
  c = zeros(n+m*d); c[:n] = 1.0
  # first try to assemble the diagonal of matrix A, and sort by its diagonal
  diagA = (X**2).sum(0) + G**2 + H**2
  isort = sorted(range(n+m*d), key=diagA.__getitem__)
  irevt = sorted(range(n+m*d), key=isort.__getitem__)
  # permute columns of X and diagonal of G and H
  X = X[:,isort]
  G = G[isort]
  H = H[isort]
  c = c[isort]
  ab = quad_prog(X, G, H, c)
  # reverse sorting permutation to get a and b and normalize
  abrevt = ab[irevt]
  a = abrevt[:n]
  b = zeros([m,d])
  for k in range(d):
    b[:,k] = abrevt[n+k*m:n+(k+1)*m]
  # compute the expeted squared residual
  finite = (ab != 0)
  Xab = dot(X[:,finite], ab[finite])**2
  Gab = (G*ab)[finite]**2
  Hab = (H*ab)[finite]**2
  er2 = Xab.sum() + Gab.sum() + Hab.sum()
  return a, b, er2



def interp_nd_grad_coef(z, x, gamma, \
                        order_set, boundary_set, boundary_zeta):
  """
  Calculate interpolation coefficients in n-dimensional spacee.
  z is the point where the interpolation is evaluated, the size of z
    is the dimension of the space.
  x is the points where the function value is available, the columns of x
    must be the same as the size of z.
  dx is the measurement error at x, its size must be same as the rows of x.
  y is the points where the function gradient is available, the columns of y
    must be the same as the size of z.
  dy is the measurement error at y, must has same shape as y.
  beta is the `magnitude' of the target function.
  gamma is the `wave number' of the target function.
    combined with beta, it provides an estimate of the derivative growth:
      f^(k) = O(beta * gamma**k)
    larger gamma = more conservative and lower order interpolation.
  order_set is the expansion order set of the Taylor expansions.
  boundary_set is the boundary of order_set.
  boundary_zeta.
  return values:
    a, b, er2 = interp_nd_coef(...)
  a and b are the interpolation coefficients.
  er2 is the expected squared residual.
  """
  # verifying arguments
  gamma = float(gamma)
  d, n, m = z.size, x.shape[0], y.shape[0]
  assert z.shape == (d,)
  assert x.shape == (n, d) and dx.shape == (n,)
  assert y.shape == (m, d) and dy.shape == (m, d) or m == 0
  assert len(order_set) > 1 and order_set[0] == (0,) * d
  # calculate interpolant at z
  d2 = ((x - z)**2).sum(1)
  if d2.min() < 1.0E-18:    # exactly matches a data point
    imatch = d2.argmin()
    if dx[imatch] == 0.0:   # and function value on that point is exact
      a = numpy.array([0]*imatch+[1]+[0]*(n-imatch-1))
      b = numpy.zeros([m, d])
      return a, b, 0.0
  # construct X = [Xa,Xb]
  N = len(order_set) - 1
  X = zeros([N, n+m*d], dtype=float)
  for i, kappa in enumerate(order_set[1:]):
    # Xa part
    X[i,:n] = gamma**sum(kappa) / factorial(kappa) * ((x - z)**kappa).prod(1)
    # Xb part
    if m > 0:
      for k in range(d):
        if kappa[k] > 0:
          kappa_p = kappa[:k] + (kappa[k] - 1,) + kappa[k+1:]
          X[i,n+k*m:n+(k+1)*m] = gamma**sum(kappa) / factorial(kappa_p) * \
                                 ((y - z)**kappa_p).prod(1)
  X *= beta
  # construct diagonal G matrix for the Lagrange residual
  G2 = zeros(n+m*d)
  for kappa, zeta in zip(boundary_set, boundary_zeta):
    Gi = gamma**sum(kappa) / factorial(kappa) * \
         ((x - z)**kappa).prod(1) * sum(zeta) / sum(kappa)
    G2[:n] += Gi**2
    if m > 0:
      for k in range(d):
        if kappa[k] > 0:
          kappa_p = kappa[:k] + (kappa[k] - 1,) + kappa[k+1:]
          zeta_p = zeta[:k] + (max(0, zeta[k] - 1),) + zeta[k+1:]
          Gi = gamma**sum(kappa) / factorial(kappa_p) * \
               ((y - z)**kappa_p).prod(1) * sum(zeta_p) / sum(kappa_p)
          G2[n+k*m:n+(k+1)*m] += Gi**2
  G = sqrt(G2) * beta
  # construct diagonal H matrix for measurement errors
  H = zeros(n+m*d)
  H[:n] = dx
  if m > 0:
    for k in range(d):
      H[n+k*m:n+(k+1)*m] = dy[:,k]
  # construct c
  c = zeros(n+m*d); c[:n] = 1.0
  # first try to assemble the diagonal of matrix A, and sort by its diagonal
  diagA = (X**2).sum(0) + G**2 + H**2
  isort = sorted(range(n+m*d), key=diagA.__getitem__)
  irevt = sorted(range(n+m*d), key=isort.__getitem__)
  # permute columns of X and diagonal of G and H
  X = X[:,isort]
  G = G[isort]
  H = H[isort]
  c = c[isort]
  ab = quad_prog(X, G, H, c)
  # reverse sorting permutation to get a and b and normalize
  abrevt = ab[irevt]
  a = abrevt[:n]
  b = zeros([m,d])
  for k in range(d):
    b[:,k] = abrevt[n+k*m:n+(k+1)*m]
  # compute the expeted squared residual
  finite = (ab != 0)
  Xab = dot(X[:,finite], ab[finite])**2
  Gab = (G*ab)[finite]**2
  Hab = (H*ab)[finite]**2
  er2 = Xab.sum() + Gab.sum() + Hab.sum()
  return a, b, er2



def calc_gamma_nd(x, fx, dfx, y, fpy, dfpy):
  """
  Estimate the `wave number' parameter gamma from data points.
  This function prints stuff.
  """
  n, m, d = x.shape[0], y.shape[0], x.shape[1]
  # calculate derivative from fx
  fpx, dfpx = zeros(n), zeros(n)
  for i in range(n):
    others = range(i) + range(i+1,n)
    d = sqrt((x - x[i,:])**2)
    fpxi = abs(fx[others,:] - fx[i,:]) / d[others]
    dfpxi = sqrt(dfx[others,:]**2 + dfx[i,:]**2) / d[others]
    imax = (fpxi / dfpxi).argmax()
    fpx[i], dfpx[i] = fpxi[imax], dfpxi[imax]
  fpy = sqrt((fpy**2).sum(1))
  dfpy = sqrt((dfpy**2).sum(1))
  # estimate beta*gamma
  ratio = ((dfpx**2).sum() + (dfpy**2).sum()) / \
          ((fpx**2).sum() + (fpy**2).sum())
  betagamma = (((fpx**2).sum() + (fpy**2).sum()) / (n+m) * exp(-ratio)) ** 0.5
  # detect discontinuity, raise gamma if needed
  if m > 0:
    dfmax = max(abs(fpx).max(), abs(fpy).max())
  else:
    dfmax = abs(fpx).max()
  betagamma = max(betagamma, 0.16 * dfmax)
  gamma = betagamma / calc_beta(fx, dfx) * 4
  print '    using gamma = ', gamma
  return gamma



def calc_res_ratio_avg_nd(beta, gamma, x, fx, dfx, y, fpy, dfpy, \
                          order_set, boundary_set, boundary_zeta):
  """
  A utility function used by calc_gamma, calculates the average ratio of real
  residual to the estimated residual, which is used to make decision in the
  bisection.
  """
  n = x.shape[0]
  n_remove = max(1, n / 5)
  i_sorted = sorted(range(n), key=lambda i:tuple(x[i]))
  var_total, var_max = 0.0, 0.0
  for i in range(0,n,n_remove):
    # interpolate each value data point using all other data points
    target = i_sorted[i:i+n_remove]
    base = i_sorted[:i] + i_sorted[i+n_remove:]
    for j in target:
      a, b, er2 = interp_nd_coef(x[j,:], x[base,:], dfx[base], y, dfpy, \
                      beta, gamma, order_set, boundary_set, boundary_zeta)
      # average the ratio of real residual to estimated residual
      resid = dot(a,fx[base]) + (b*fpy).sum() - fx[j]
      var = resid**2 / (er2 + dfx[j]**2)
      var_total += var
      var_max = max(var_max, var)
  return var_total / n



def calc_gamma_bisect_nd(x, fx, dfx, y, fpy, dfpy, \
                  order_set, boundary_set, boundary_zeta):
  """
  Estimate the `wave number' parameter gamma from data points.
  This function prints stuff.
  """
  beta = calc_beta(fx, dfx)
  print '    using beta = ', beta
  # upper and lower bounds
  delta_min, delta_max = numpy.inf, 0.0
  for xyi in list(x) + list(y):
    d2x = ((xyi - x)**2).sum(1)
    if y.size > 0:
      d2y = ((xyi - y)**2).sum(1)
      delta_max = max(delta_max, max(d2x.max(), d2y.max()))
    else:
      delta_max = max(delta_max, d2x.max())
  for i, xi in enumerate(x):
    if i > 0:
      d2x = ((xi - x[:i,:])**2).sum(1)
      delta_min = min(delta_min, d2x.min())
  for i, yi in enumerate(y):
    if i > 0:
      d2y = ((yi - y[:i,:])**2).sum(1)
      delta_min = min(delta_min, d2y.min())
  delta_min, delta_max = sqrt(delta_min), sqrt(delta_max)
  assert delta_max > delta_min
  gamma_min = 1. / delta_max
  gamma_max = pi / delta_min
  # logorithmic bisection for gamma
  while gamma_max / gamma_min > 1.1:
    print '    bisecting [', gamma_min, ',', gamma_max, '] for gamma...'
    gamma_mid = sqrt(gamma_max * gamma_min)
    res_ratio = calc_res_ratio_avg_nd(beta, gamma_mid, x, fx, dfx, \
                y, fpy, dfpy, order_set, boundary_set, boundary_zeta)
    if res_ratio < 1.0:
      gamma_max = gamma_mid
    else:
      gamma_min = gamma_mid
  # final selected gamma
  gamma_mid = sqrt(gamma_max * gamma_min)
  print '    using gamma = ', gamma_mid
  return gamma_mid



def interp_nd(z, x, fx, dfx=None, y=None, fpy=None, dfpy=None, \
              compute_dfz=False, order_set=None):
  """
  Wang interpolation in d-dimensional space.
  z is the points where the interpolation is evaluated.
  x is the points where the function value is available and given by fx, dfx
    estimates the error in the function values, default to 0 (exact value).
  y is the points where the function gradient is available and given by fpy,
    dfpy estimates the error in the gradient values, default to 0
    (exact gradient).
  compute_dfz indicates whether an estimate of error in the interpolation
    approximation is returned.  If False (default), "fz = interp_1d(...)";
    if true, "fz, dfz = interp_1d(...)"
  """
  # verifying and handling arguments
  d, n = z.shape[1], x.shape[0]
  assert x.shape == (n,d) and fx.shape == (n,)
  if dfx is None:
    dfx = zeros(fx.shape)
  else:
    assert dfx.shape == (n,)
  if y is None:
    assert fpy is None and dfpy is None
    m = 0; y = zeros([0,d]); fpy = zeros([0,d]); dfpy = zeros([0,d])
  else:
    m = y.shape[0]
    assert fpy is not None and fpy.shape == (m,d)
    if dfpy is None:
      dfpy = zeros(fpy.shape)
    else:
      assert dfpy.shape == (m,d)
  # determine expansion order set, its boundary set and zeta
  if order_set is None:
    k = 0; order_set = []
    while len(order_set) < min(n + m*d, 100):
      k += 1
      order_set = max_order_set(d, k)
  else:
    assert len(order_set) > 1 and order_set[0] == (0,)*d
  boundary_set, boundary_zeta = boundary_set_and_zeta(order_set)
  # calculate beta and gamma
  t0 = time.time()
  beta = calc_beta(fx, dfx)
  gamma = calc_gamma_nd(x, fx, dfx, y, fpy, dfpy, \
                        order_set, boundary_set, boundary_zeta)
  print 'time: ', time.time() - t0
  # interpolation for each z[i]
  t0 = time.time()
  fz, dz = [], []
  for zi in z:
    a, b, er2 = interp_nd_coef(zi, x, dfx, y, dfpy, beta, gamma, \
                               order_set, boundary_set, boundary_zeta)
    fz.append(dot(a, fx) + (b*fpy).sum())
    dz.append(sqrt(er2))
  print 'time: ', time.time() - t0
  if compute_dfz:
    return numpy.array(fz), numpy.array(dz)
  else:
    return numpy.array(fz)



