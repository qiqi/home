"""
Qiqi Wang  January 2009  Interpolation and regression in one-dimensional domain

reference:
  Q Wang et al. A Rational interpolation Scheme with Super-Polynomial Rate of
    Convergence.  Submitted to SIAM Journal of Numerical Analysis.
"""


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
  For a list of integer, factorial([a,b,c]) = a! * b! * c!
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
  Solve a lower triangular system
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
  Solve an upper triangular system
  """
  n = b.size
  assert R.shape == (n,n)
  x = zeros(n, dtype=R.dtype)
  for i in range(n-1,-1,-1):
    x[i] = (b[i] - dot(x[i+1:], R[i,i+1:])) / R[i,i]
    if not numpy.isfinite(x[i]):
      x[i] = 0.0
  return x



def quad_prog_wang(X, E, c):
  """
  Solve the constraint quadratic programming
    min a^T A a, where A = X^T X + E^2, s.t. a^T c = 1
  This only works when c is 1-D.  Use quad_prog_golub for multiple constraints.
  To get best performance, order columns of X and E such that the diagonal
    of A is of increasing size.
  """
  # filter out infinite columns and rows
  D = (X**2).sum(0) + E**2
  finite = numpy.isfinite(D)
  X = X[:,finite]
  E = E[finite]
  # QR decomposition of X, R is now cholesky factor of A
  N, n = X.shape
  A = zeros([N+n, n])
  A[:N,:] = X
  A[range(N,N+n), range(n)] = E
  R = numpy.linalg.qr(A, mode='r')
  # solve...
  ab = zeros(c.size, dtype=X.dtype)
  tmp = solve_L(R.transpose(), c[finite])
  ab[finite] = solve_R(R, tmp)
  ab /= dot(c, ab)
  return asarray(ab)



def quad_prog_golub(X, E, C):
  """
  Solve the constraint quadratic programming with
    min a^T A a, where A = X^T X + E^2, s.t. a^T C = 1
  Algorithm in Gene Golub's book.
  To get best performance, order columns of X and E such that the diagonal
    of A is of increasing size.
  """
  # Householder transformation for constraints C
  if C.ndim == 1:
    C = C.reshape([1, C.size])
  assert X.shape[1] == E.shape[0] and X.ndim == 2 and E.ndim == 1
  assert C.ndim == 2 and C.shape[1] == E.size and C.shape[0] < C.shape[1]
  l = C.shape[0]
  # filter out infinite columns and rows
  D = (X**2).sum(0) + E**2
  finite = numpy.isfinite(D)
  X = X[:,finite]
  E = E[finite]
  C = C[:,finite]
  # prepare householder
  v = zeros(C.shape)
  E = diag(E)
  for ic in range(l):
    u = C[ic,ic:].copy()
    u[0] += sign(C[ic,ic]) * sqrt((C[ic,ic:]**2).sum())
    v[ic,ic:] = u / sqrt((u**2).sum())
    # for C
    cv = dot(C[:,ic:],v[ic,ic:])
    for i in range(ic,l):
      C[i,ic:] -= 2 * cv[i] * v[ic,ic:]
    # for X
    xv = dot(X[:,ic:], v[ic,ic:])
    for i in range(xv.size):
      X[i,ic:] -= 2 * xv[i] * v[ic,ic:]
    # householder for E
    hgv = dot(E[:,ic:], v[ic,ic:])
    for i in range(hgv.size):
      E[i,ic:] -= 2 * hgv[i] * v[ic,ic:]
  # QR decomposition
  N, n = X.shape
  A = zeros([N+n, n])
  A[:N,:] = X[:,:]
  A[N:, :] = E[:,:]
  Q, R = numpy.linalg.qr(A[:,l:], mode='full')
  # solve for constraints
  at = zeros(n)
  rhs1 = zeros(l); rhs1[0] = 1.0
  at[:l] = linalg.solve(C[:,:l],rhs1)
  # solve for least squares
  rhs2 = -dot(Q.transpose(), dot(A[:,:l], at[:l]))
  at[l:] = solve_R(R, rhs2)
  # invert the householder
  for ic in range(l-1, -1, -1):
    at[ic:] -= 2 * v[ic,ic:] * dot(at, v[ic,:])
  a = zeros(finite.size)
  a[finite] = at
  return asarray(a)



class interp_1d:
  """
  Interpolation in 1D.
  xv is the ``value data points'', i.e. the points where the function value
     is available and given by fxv; dfxv estimates the standard deviation of
     the error in the function values; the default of dfxv is 0 (fxv is exact).
  xg is the ``gradient data points'', i.e. points where the function gradient
     is available and given by fpxg; dfpxg estimates the standard devisiton of
     the error in the gradient values; the default of dfpxg is 0
     (dfxg is exact).
  beta is the `magnitude' of the target function, can be automatically
     calculated.
  gamma is the `wave number' of the target function, can be automatically
     calculated.
     Combined with beta, it provides an estimate of the derivative growth:
     f^(k) = O(beta * gamma**k)
     Larger gamma yields more conservative, more robust and lower order
     interpolation.
  N is the order of the Taylor expansion, can be automatically calculated.
     Smaller N yields lower order interpolation.  Numerical instability may
     occur when N is too large.
  l is the polynomial order.  The interpolant is forced to interpolate order
     l-1 polynomials exactly.  l=1 is the most robust, higher l makes a
     difference only when gamma is large, or when data is sparse and oscilatory
     if gamma is automatically calculated.
  verbose is the verbosity level.  0 is silent.
  Reference:
  * Q.Wang et al. A Rational Interpolation Scheme with Super-polynomial Rate
    of Convergence.
  """

  def __init__(self, xv, fxv, dfxv=None, xg=None, fpxg=None, dfpxg=None, \
               beta=None, gamma=None, N=None, l=1, verbose=1):
    """
    __init__(self, xv, fxv, dfxv=None, xg=None, fpxg=None, dfpxg=None,
             beta=None, gamma=None, N=None, l=1)
    Instantiation function, see class documentation for arguments.
    fxv must has same size as xv.
    dfxv must has same size as xv, or None for default (all 0).
    fpxg must be None if xg is None, or has same size as xg if xg is not None.
    dfpxg must be None if xg is None; if xg is not None it must has same size
      as xg, or None for default (all 0).
    When beta and gamma must be both None or both given.  When they are
      None, their values are automatically calculated.  The calculation of
      gamma may take a long time if the number of datapoints is large.
    """
    
    assert verbose == 0 or verbose == 1 or verbose == 2
    self.verbose = verbose

    # verify and save value data points
    assert xv.ndim == fxv.ndim == 1
    assert xv.size == fxv.size
    self.xv = copy.copy(xv)
    self.fxv = copy.copy(fxv)
    if dfxv is None:
      self.dfxv = zeros(xv.shape)
    else:
      assert dfxv.shape == xv.shape
      self.dfxv = copy.copy(dfxv)

    # verify and save gradient data points
    if xg is None:
      assert fpxg is None and dfpxg is None
      self.xg = zeros(0)
      self.fpxg = zeros(0)
      self.dfpxg = zeros(0)
    else:
      assert fpxg is not None and fpxg.shape == xg.shape
      self.xg = copy.copy(xg)
      self.fpxg = copy.copy(fpxg)
      if dfpxg is None:
        self.dfpxg = zeros(xg.shape)
      else:
        assert dfpxg.shape == xg.shape
        self.dfpxg = copy.copy(dfpxg)

    # check and automatically calculate N
    self.nv = self.xv.size
    self.ng = self.xg.size
    self.n = self.nv + self.ng
    if N is None:
      self.N = min(self.n, 200)
    else:
      self.N = N

    assert int(l) == l
    self.l = int(l)

    # automatically calculate beta and gamma
    if beta is None:
      assert gamma is None
      self.beta = self.calc_beta()
      self.gamma = self.calc_gamma()
    else:
      self.beta = float(beta)
      self.gamma = float(gamma)



  def interp_matrices(self, x):
    """
    Calculate the matrices X, E and C for interpolation scheme
    x is the point where the interpolation is evaluated.
    """
    assert isinstance(x, float)
    N, n, l, nv, ng = self.N, self.n, self.l, self.nv, self.ng

    # construct X = [Xv,Xg]
    X = zeros([N, n], dtype=float)
    for i in range(N):
      X[i,:nv] = self.gamma**(i+1) / factorial(i+1) * (self.xv - x) ** (i+1)
      X[i,nv:] = self.gamma**(i+1) / factorial(i) * (self.xg - x) ** i
    X *= self.beta

    # construct diagonal Er matrix for the Lagrange residual
    Er = zeros(n)
    Er[:nv] = self.gamma**(N+1) / factorial(N+1) * (self.xv - x)**(N+1)
    Er[nv:] = self.gamma**(N+1) / factorial(N) * (self.xg - x)**N
    Er *= self.beta
    # construct diagonal Ee matrix for measurement errors
    Ee = zeros(n)
    Ee[:nv] = self.dfxv
    Ee[nv:] = self.dfpxg
    # construct E
    E = sqrt(Er**2 + Ee**2)

    # construct C
    C = zeros([l,n])
    C[0,:nv] = 1.0
    C[0,nv:] = 0.0
    for i in range(1, l):
      C[i,:nv] = (self.xv-x)**i
      C[i,nv:] = i * (self.xg-x)**(i-1)
    return X, E, C



  def interp_coef(self, x):
    """
    Calculate interpolation coefficients in 1D.
    x is the point where the interpolation is evaluated.
    return values:
      av, ag, er2 = obj.interp_coef(x)
    av and ag are the interpolation coefficients for value data points and
      gradient data points, respectively
    er2 is the estimated variance of the interpolation residual.
    """
    x = float(x)
    # calculate interpolant at z
    if abs(self.xv - x).min() < 1.0E-12:   # exactly matches a data point
      imatch = abs(self.xv - x).argmin()
      if self.dfxv[imatch] == 0.0:   # and function value on that point is exact
        av = numpy.array([0]*imatch+[1]+[0]*(self.nv-imatch-1))
        ag = numpy.zeros(self.ng)
        return av, ag, 0.0
    # construct X = [Xv,Xg]
    X, E, C = self.interp_matrices(x)
    # first try to assemble the diagonal of matrix A, and sort by its diagonal
    diagA = (X**2).sum(0) + E**2
    isort = sorted(range(self.n), key=diagA.__getitem__)
    irevt = sorted(range(self.n), key=isort.__getitem__)
    # permute columns of X and diagonal of G and H
    X = X[:,isort]
    E = E[isort]
    C = C[:,isort]
    if C.shape[0] == 1:
      a = quad_prog_wang(X, E, C[0,:])
    else:
      a = quad_prog_golub(X, E, C)
    # reverse sorting permutation to get a and b and normalize
    arevt = a[irevt]
    av = arevt[:self.nv]
    ag = arevt[self.nv:]
    # compute the expeted squared residual
    finite = (a != 0)
    Xa = dot(X[:,finite], a[finite])**2
    Ea = (E*a)[finite]**2
    er2 = Xa.sum() + Ea.sum()
    return av, ag, er2



  def calc_beta(self):
    """
    Estimate the `magnitude' parameter beta from data points.
    """
    assert self.fxv.ndim == 1 and self.fxv.shape == self.dfxv.shape
    f_bar = self.fxv.mean()
    ratio = (self.dfxv**2).sum() / ((self.fxv - f_bar)**2).sum() * \
            float(self.nv-1) / float(self.nv)
    beta = sqrt(((self.fxv - f_bar)**2).sum() / (self.nv-1) * exp(-ratio))
    return beta



  def calc_res_ratio_avg(self, beta, gamma):
    """
    A utility function used by calc_gamma, calculates the average ratio of real
    residual to the estimated residual, which is used to make decision in the
    bisection.
    """
    n_remove = max(1, self.nv / 5)
    i_sorted = sorted(range(self.nv), key=self.xv.__getitem__)
    var_total = 0.0
    for i in range(0,self.nv,n_remove):
      # interpolate each value data point using all other data points
      target = i_sorted[i:i+n_remove]
      base = i_sorted[:i] + i_sorted[i+n_remove:]
      subinterp = interp_1d(self.xv[base], self.fxv[base], self.dfxv[base], \
                            self.xg, self.fpxg, self.dfpxg, beta, gamma, \
                            self.N, self.l, self.verbose)
      for j in target:
        av, ag, er2 = subinterp.interp_coef(self.xv[j])
        # average the ratio of real residual to estimated residual
        resid = dot(av,self.fxv[base]) + dot(ag,self.fpxg) - self.fxv[j]
        var_total += resid**2 / (er2 + self.dfxv[j]**2)
    return var_total / self.nv



  def calc_gamma(self):
    """
    Estimate the `wave number' parameter gamma from data points.
    This function prints stuff when self.verbose > 0.
    """
    assert isinstance(self.beta, float)
    # upper and lower bounds
    sorted_xv = numpy.array(sorted(self.xv))
    sorted_xg = numpy.array(sorted(self.xg))
    if self.xg.size == 0:
      delta_max = self.xv.max() - self.xv.min()
      delta_min = (sorted_xv[1:] - sorted_xv[:-1]).min()
    else:
      delta_max = max(self.xv.max(), self.xg.max()) - \
                  min(self.xv.min(), self.xg.min())
      delta_min = min((sorted_xv[1:] - sorted_xv[:-1]).min(), \
                      (sorted_xg[1:] - sorted_xg[:-1]).min())
    assert delta_max > delta_min
    gamma_min = 1. / delta_max
    gamma_max = pi / delta_min
    # logorithmic bisection for gamma
    while gamma_max / gamma_min > 1.1:
      if self.verbose > 0:
        print '    bisecting [', gamma_min, ',', gamma_max, '] for gamma...'
      gamma_mid = sqrt(gamma_max * gamma_min)
      res_ratio = self.calc_res_ratio_avg(self.beta, gamma_mid)
      if res_ratio < 1.0:
        gamma_max = gamma_mid
      else:
        gamma_min = gamma_mid
    # final selected gamma
    gamma_mid = sqrt(gamma_max * gamma_min)
    if self.verbose > 0:
      print '    using gamma = ', gamma_mid
    return gamma_mid



  def interp(self, x, compute_df=False):
    """
    Interpolation in 1D.
    x is the point (a single number) or points (a list of numbers) where the
      interpolation is evaluated.
    compute_df indicates whether an estimated standard deviation of the error
      in the interpolation approximation is also returned.
    Usage:
      fx = interp(x, compute_df=False)";
      fx, df = interp(x, compute_df=True)"
    """
    # evaluate interpolant value at a single point
    if isinstance(x, (float, int)):
      av, ag, er2 = self.interp_coef(x)
      fx = dot(av, self.fxv) + dot(ag, self.fpxg)
      dfx = sqrt(er2)
      if compute_df:
        return fx, dfx
      else:
        return fx
    # at multiple points
    else:
      fx, dfx = [], []
      for xi in x:
        av, ag, er2 = self.interp_coef(xi)
        fx.append(dot(av, self.fxv) + dot(ag, self.fpxg))
        dfx.append(sqrt(er2))
      if compute_df:
        return numpy.asarray(fx), numpy.asarray(dfx)
      else:
        return numpy.asarray(fx)

