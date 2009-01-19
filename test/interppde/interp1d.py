import copy
import sys
import pickle
import time

import numpy
import pylab
from numpy import zeros, ones, eye, kron, linalg, dot, exp, sqrt, diag, pi, \
                  asarray, sign



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
  # filter out infinite columns and rows
  D = (X**2).sum(0) + G**2 + H**2
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
  return asarray(ab)



def quad_prog_new(X, G, H, c):
  """
  Solve the constraint quadratic programming with
  min ab^T A ab, where A = X^T X + G^2 + H^2, s.t. ab^T c = 1
  """
  # Householder transformation for constraints c
  if c.ndim == 1:
    c = c.reshape([1, c.size])
  assert c.ndim == 2 and c.shape[1] == H.size and c.shape[0] < c.shape[1]
  nc = c.shape[0]
  # filter out infinite columns and rows
  D = (X**2).sum(0) + G**2 + H**2
  finite = numpy.isfinite(D)
  X = X[:,finite]
  H = H[finite]
  G = G[finite]
  c = c[:,finite]
  # prepare householder
  v = zeros(c.shape)
  HG = diag(sqrt(H**2 + G**2))
  for ic in range(nc):
    u = c[ic,ic:].copy()
    u[0] += sign(c[ic,ic]) * sqrt((c[ic,ic:]**2).sum())
    v[ic,ic:] = u / sqrt((u**2).sum())
    # for c
    cv = dot(c[:,ic:],v[ic,ic:])
    for i in range(ic,nc):
      c[i,ic:] -= 2 * cv[i] * v[ic,ic:]
    # for X
    xv = dot(X[:,ic:], v[ic,ic:])
    for i in range(xv.size):
      X[i,ic:] -= 2 * xv[i] * v[ic,ic:]
    # householder for HG
    hgv = dot(HG[:,ic:], v[ic,ic:])
    for i in range(xv.size):
      HG[i,ic:] -= 2 * hgv[i] * v[ic,ic:]
  # QR decomposition
  n1, n2 = X.shape
  A = zeros([n1+n2, n2])
  A[:n1,:] = X[:,:]
  A[n1:, :] = HG[:,:]
  Q, R = numpy.linalg.qr(A[:,nc:], mode='full')
  # solve for constraints
  abt = zeros(n2)
  rhs1 = zeros(nc); rhs1[0] = 1.0
  abt[:nc] = linalg.solve(c[:,:nc],rhs1)
  # solve for least squares
  rhs2 = -dot(Q.transpose(), dot(A[:,:nc], abt[:nc]))
  abt[nc:] = solve_R(R, rhs2)
  # invert the householder
  for ic in range(nc-1, -1, -1):
    abt[ic:] -= 2 * v[ic,ic:] * dot(abt, v[ic,:])
  ab = zeros(finite.size)
  ab[finite] = abt
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
  c = zeros([2,n+m])
  c[0,:n] = 1.0; c[0,n:] = 0.0
  c[1,:n] = x-z; c[1,n:] = 1.0
  # first try to assemble the diagonal of matrix A, and sort by its diagonal
  diagA = (X**2).sum(0) + G**2 + H**2
  isort = sorted(range(n+m), key=diagA.__getitem__)
  irevt = sorted(range(n+m), key=isort.__getitem__)
  # permute columns of X and diagonal of G and H
  X = X[:,isort]
  G = G[isort]
  H = H[isort]
  c = c[:,isort]
  ab = quad_prog_new(X, G, H, c)
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
  x[ix] is the point where the interpolant gradient is evaluated.
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
      a, b, er2 = interp_1d_coef(x[j], x[base], dfx[base], \
                                      y, dfpy, beta, gamma)
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
  Interpolation in 1D.
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

