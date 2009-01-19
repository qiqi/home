class Lorenz (object):
   def __init__(self, sigma, rho, beta):
      self._sigma = sigma
      self._rho = rho
      self._beta = beta

   def step(self, x, y, z, dt):
      dx = dt * self._sigma * (y - x)
      dy = dt * (self._rho * x - x * z - y)
      dz = dt * (x * y - self._beta * z)
      return x + dx, y + dy, z + dz

   def diff(self, x, y, z, dx, dy, dz, dt):
      ddx = dt * self._sigma * (dy - dx)
      ddy = dt * (self._rho * dx - x * dz - z * dx - dy)
      ddz = dt * (x * dy + y * dx - self._beta * dz)
      return dx + ddx, dy + ddy, dz + ddz

   def adj(self, x, y, z, xadj, yadj, zadj, dt):
      dxadj = dt * (- self._sigma * xadj + (self._rho - z) * yadj + y * zadj)
      dyadj = dt * (self._sigma * xadj - yadj + x * zadj)
      dzadj = dt * (- x * yadj - self._beta * zadj)
      return xadj + dxadj, yadj + dyadj, zadj + dzadj
