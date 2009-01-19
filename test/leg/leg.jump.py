import logging

import numpy
from numpy import eye, diag, zeros, asmatrix, asarray, sin, cos, dot, squeeze
from numpy.linalg import inv, solve, norm
import pylab

logging.basicConfig()
_logger = logging.Logger('LEG')

class MonoLeg:
   # Simulate a robot with a body and a single leg consisting of two pieces.
   # Totally three pieces, 5 degrees of freedom in 2D:
   #    one position in 2D (2) and three angles (3).
   # We simulate it with 9 degrees of freedom:
   #    three positions in 2D (6) and 3 angles (3).
   # and pose 4 constraints:
   #    position of two joints (4).
   # The lagrange multiplier of the 4 constraints are the (x,y) forces at the
   # two joint, denoted as lambdas.
   # The controls in the system is the bending velocities.

   def __init__(self, mass, length, inertia, g_std=9.8, n_newton_iter=5):
      # Construct linear system
      #    M dV/dt = G + F [lambdas]
      #    dX/dt = V
      #    Constraint(X,V,controls) == 0
      assert len(mass) == len(length) == len(inertia) == 3
      self._mass = asarray(mass, dtype=float)
      self._length = asarray(length, dtype=float)
      self._inertia = asarray(inertia, dtype=float)
      self._g_std = g_std
      self._n_newton_iter = n_newton_iter
      # state variables: V and X
      self._V0 = numpy.zeros(9)
      self._V = numpy.zeros(9)
      self._X0 = numpy.zeros(9)
      self._X = numpy.zeros(9)
      # parts of X, just to be handy
      self._x_pos = self._X[0:9:3]
      self._y_pos = self._X[1:9:3]
      self._a_pos = self._X[2:9:3]
      self._x_vel = self._V[0:9:3]
      self._y_vel = self._V[1:9:3]
      self._a_vel = self._V[2:9:3]
      # initialize, standing upright
      self._x_pos[0:3] = [0.0, 0.0, 0.0]
      self._y_pos[0:3] = [0.5 * length[0] + length[1] + length[2], \
                          0.5 * length[1] + length[2], 0.5 * length[2]]
      self._a_pos[0:3] = [0.0, 0.0, 0.0]
      # Calculate M, G
      self._M = diag([mass[0], mass[0], inertia[0], \
                      mass[1], mass[1], inertia[1], \
                      mass[2], mass[2], inertia[2]])
      self._M_inv = inv(self._M)
      self._G = -self._g_std * \
                numpy.array([0,mass[0],0, 0,mass[1],0, 0,mass[2],0])
      # initialize F matrix
      self._F = zeros((9,8))
      self._F[0:3,0:3], self._F[3:6,0:3] = eye(3), -eye(3)
      self._F[3:6,3:6], self._F[6:9,3:6] = eye(3), -eye(3)
      self._F[6:8,6:8] = eye(2)
      #self.calc_F_matrix()
      # initialize C matrix
      self._C = zeros((8,9))
      self._C[0:2,0:2], self._C[0:2,3:5] = eye(2), -eye(2)
      self._C[2:4,3:5], self._C[2:4,6:8] = eye(2), -eye(2)
      self._C[4,2], self._C[4,5], self._C[5,5], self._C[5,8] = 1,-1,1,-1
      self._C[6:8,6:8] = eye(2)
      # initialize lambda and control
      self._lambda = zeros(8)
      self._control = zeros(2)

   def step(self, dt):
      self.calc_F_matrix()
      # save old velocity and update new one (predictor)
      self._V0[0:9] = self._V
      self._V += dt * dot(self._M_inv, self._G)
      # update new position
      self._X0[0:9] = self._X
      self._X += dt * (self._V0 + self._V) / 2.0
      # calculate lagrange multiplier to correct joints (corrector)
      self.corrector(dt)

   def corrector(self, dt):
      r = zeros(8)
      self._lambda[0:8] = 0.0
      for i_iter in range(self._n_newton_iter):
         # calculate residual
         cos_arms = 0.5 * self._length * cos(self._a_pos)
         sin_arms = 0.5 * self._length * sin(self._a_pos)
         r[0] = self._x_pos[0] - self._x_pos[1] - sin_arms[0] - sin_arms[1]
         r[1] = self._y_pos[0] - self._y_pos[1] - cos_arms[0] - cos_arms[1]
         r[2] = self._x_pos[1] - self._x_pos[2] - sin_arms[1] - sin_arms[2]
         r[3] = self._y_pos[1] - self._y_pos[2] - cos_arms[1] - cos_arms[2]
         r[4] = self._a_vel[0] - self._a_vel[1] - self._control[0]
         r[5] = self._a_vel[1] - self._a_vel[2] - self._control[1]
         r[0:4] /= 0.5 * dt**2
         r[4:6] /= dt
         self.calc_C_matrix()
         J = asmatrix(self._C[0:6,:]) * asmatrix(self._M_inv) * \
             asmatrix(self._F[:,0:6])
         self._lambda[0:6] -= solve(J, r[0:6])
         _logger.debug('Air Newton iter %d norm = %f' % (i_iter, norm(r[0:6])))
         self._V[0:9] = self._V0 + dt * dot(self._M_inv, self._G + \
                        dot(self._F[:,0:6], self._lambda[0:6]))
         self._X[0:9] = self._X0 + dt * (self._V0 + self._V) / 2.0
      # determine if foot is in air or on ground
      if self._y_pos[2] - cos_arms[2] < 0:
         # calculate foot x position
         x1, y1 = self._x_pos[2] - sin_arms[2], self._y_pos[2] - cos_arms[2]
         x0 = self._X0[6] - 0.5 * self._length[2] * sin(self._X0[8])
         y0 = self._X0[7] - 0.5 * self._length[2] * cos(self._X0[8])
         if y1 - y0 != 0:
             foot_pos = x0 - (x1 - x0) / (y1 - y0) * y0
         else:
             foot_pos = x0
         for i_iter in range(self._n_newton_iter):
            cos_arms = 0.5 * self._length * cos(self._a_pos)
            sin_arms = 0.5 * self._length * sin(self._a_pos)
            r[0] = self._x_pos[0] - self._x_pos[1] - sin_arms[0] - sin_arms[1]
            r[1] = self._y_pos[0] - self._y_pos[1] - cos_arms[0] - cos_arms[1]
            r[2] = self._x_pos[1] - self._x_pos[2] - sin_arms[1] - sin_arms[2]
            r[3] = self._y_pos[1] - self._y_pos[2] - cos_arms[1] - cos_arms[2]
            r[4] = self._a_vel[0] - self._a_vel[1] - self._control[0]
            r[5] = self._a_vel[1] - self._a_vel[2] - self._control[1]
            r[6] = self._x_pos[2] - sin_arms[2] - foot_pos
            r[7] = self._y_pos[2] - cos_arms[2]
            r[0:4] /= 0.5 * dt**2
            r[6:8] /= 0.5 * dt**2
            r[4:6] /= dt
            _logger.debug('Land Newton iter %d norm = %f' % \
                          (i_iter, norm(r[0:8])))
            self.calc_C_matrix()
            J = asmatrix(self._C) * asmatrix(self._M_inv) * asmatrix(self._F)
            self._lambda -= solve(J, r[0:8])
            self._V[0:9] = self._V0 + dt * dot(self._M_inv, self._G + \
                           dot(self._F, self._lambda))
            self._X[0:9] = self._X0 + dt * (self._V0 + self._V) / 2.0
      
   def calc_C_matrix(self):
      # calculate the constraint Jacobian matrix C, used in Newton iteration
      cos_arms = 0.5 * self._length * cos(self._a_pos)
      sin_arms = 0.5 * self._length * sin(self._a_pos)
      self._C[0:2,2] = -cos_arms[0], sin_arms[0]
      self._C[0:2,5] = -cos_arms[1], sin_arms[1]
      self._C[2:4,5] = -cos_arms[1], sin_arms[1]
      self._C[2:4,8] = -cos_arms[2], sin_arms[2]
      self._C[6:8,8] = -cos_arms[2], sin_arms[2]

   def calc_F_matrix(self):
      # calculate the matrix before the lagrange multipliers
      # (the forces and angular moments at the two joints)
      cos_arms = 0.5 * self._length * cos(self._a_pos)
      sin_arms = 0.5 * self._length * sin(self._a_pos)
      self._F[2,0] = -cos_arms[0]
      self._F[2,1] = sin_arms[0]
      self._F[5,0] = cos_arms[1]
      self._F[5,1] = -sin_arms[1]
      self._F[5,3] = -cos_arms[1]
      self._F[5,4] = sin_arms[1]
      self._F[8,3] = cos_arms[2]
      self._F[8,4] = -sin_arms[2]
      self._F[8,6] = -cos_arms[2]
      self._F[8,7] = sin_arms[2]

   def plot(self, xlim=None, ylim=None):
      pylab.cla()
      x, y = zeros((2,3)), zeros((2,3))
      x[0,0:3] = self._x_pos + 0.5 * self._length * sin(self._a_pos)
      y[0,0:3] = self._y_pos + 0.5 * self._length * cos(self._a_pos)
      x[1,0:3] = self._x_pos - 0.5 * self._length * sin(self._a_pos)
      y[1,0:3] = self._y_pos - 0.5 * self._length * cos(self._a_pos)
      pylab.plot(x, y, '-+')
      pylab.axis('equal')
      if xlim is not None:
          pylab.plot(xlim, [0,0], ':k')
          pylab.xlim(xlim)
      if ylim is not None:
          pylab.ylim(ylim)

