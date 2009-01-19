import numpy
from numpy import zeros, ones, linalg, diag, isfinite

def lagrange_1d(y, x, fx):
  """
  Lagrange polynomial interpolation in 1D.
  """
  fy = zeros(y.shape, dtype=float)
  for xi, fxi in zip(x, fx):
    fyi = ones(y.shape, dtype=float) * fxi
    for xj in x:
      if xi != xj:
        fyi *= (y - xj) / (xi - xj)
    fy += fyi
  return fy



def linear_1d(y, x, fx):
  """
  Linear interpolation in 1D
  """
  fy = []
  for yi in y:
    xp = x[x >= yi].min()
    xm = x[x <= yi].max()
    fp = fx[x == xp][0]
    fm = fx[x == xm][0]
    if xp - xm <= 0:
      fy.append(fp)
    else:
      fy.append((fp * (yi - xm) + fm * (xp - yi)) / (xp - xm))
  return numpy.array(fy)



def spline_1d(y, x, fx):
  """
  Linear interpolation in 1D
  """
  n = x.size
  isort = sorted(range(n), key=x.__getitem__)
  x = x[isort]; fx = fx[isort]
  h = x[1:] - x[:-1]
  # solve for spline coefficients
  A = diag(2 * (h[1:] + h[:-1])) + diag(h[1:-1], 1) + diag(h[1:-1], -1)
  rhs = 6 * ((fx[2:] - fx[1:-1]) / h[1:] - (fx[1:-1] - fx[:-2]) / h[:-1])
  z = zeros(n)
  z[1:-1] = linalg.solve(A, rhs)
  # calculate spline
  fy = []
  for yi in y:
    if abs(x - yi).min() < 1.0E-12:
      fy.append(fx[abs(x - yi).argmin()])
    else:
      im = (x < yi).sum() - 1
      if im < 0: im = 0
      if im > n-2: im = n-2
      ip = im + 1
      dxp = x[ip] - yi
      dxm = yi - x[im]
      hi = h[im]
      fy.append((z[ip] * dxm**3 + z[im] * dxp**3) / (6*hi) + \
                (fx[ip]/hi - hi/6*z[ip]) * dxm + (fx[im]/hi - hi/6*z[im]) * dxp)
  return numpy.array(fy)


def floater_hormann(y, x, fx, d):
  """
  Floater-Hormann rational interpolation with no poles
  """
  assert x.ndim == y.ndim == 1
  assert x.shape == fx.shape
  n = x.size
  assert 1 <= d <= n
  # sort x
  isort = sorted(range(n), key=x.__getitem__)
  x, fx = x[isort], fx[isort]
  # construct polynomial interpolations
  P = zeros([n-d+1, y.size])
  for i in range(n-d+1):
    P[i,:] = lagrange_1d(y, x[i:i+d], fx[i:i+d])
  # construct lambda
  Lambda = zeros([n-d+1, y.size])
  for i in range(n-d+1):
    Lambda[i,:] = (-1)**i
    for j in range(i, i+d):
      Lambda[i,y!=x[j]] /= (y-x[j])[y!=x[j]]
      Lambda[i,y==x[j]] = numpy.inf
  # modify infinities
  for i in range(y.size):
    if not isfinite(Lambda[:,i]).all():
      jinf = isfinite(Lambda[:,i]).argmin()
      Lambda[:,i] = 0.0
      Lambda[jinf,i] = 1.0
  return (Lambda * P).sum(0) / Lambda.sum(0)



def floater_hormann_adaptive(y, x, fx):
  """
  Floater-Hormann rational interpolation with adaptive d
  """
  assert x.ndim == y.ndim == 1
  assert x.shape == fx.shape
  n = x.size
  # sort x
  isort = sorted(range(n), key=x.__getitem__)
  x, fx = x[isort], fx[isort]
  # find the right order
  linf = zeros(n-2)
  for d in range(1, n-1):
    for i in range(n-1):
      xx = numpy.array(list(x[:i]) + list(x[i+2:]))
      fxx = numpy.array(list(fx[:i]) + list(fx[i+2:]))
      fxi = floater_hormann(x[i:i+2], xx, fxx, d)
      err = max(abs(fxi[0] - fx[i]), abs(fxi[1] - fx[i+1]))
      linf[d-1] = max(linf[d-1], err)
    print '   trying d = ', d, ', L_inf error = ', linf[d-1]
    if d >= 3 and linf[d-1] > linf[d-2] > linf[d-3]:
      linf[d:] = numpy.inf
      break
  d = linf.argmin() + 1
  print '   using d = ', d
  return floater_hormann(y, x, fx, d)

