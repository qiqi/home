import sys

import numpy
from numpy import zeros, logical_and

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from continent import continent, slope
from weather import temperature, precipitation
from vegetation import vegetation

SCALE = 10.0

def initScene():
   print 'Loading...'
   # heightmap
   h = continent(4,5,0.465,8)    # penesulas
   # h = continent(7,7,0.0,8)      # mountainous
   # h = continent(13,13,0.0,8)    # good plains
   # h = continent(24,65,0.0,8)    # sahara
   h *= SCALE

   summer, winter = temperature(h, 10, 60)
   precip = precipitation(h, 10, 60)
   dh = slope(h)
   r, g, b = vegetation(h, dh, summer, winter, precip)

   # h[h == 0] = - 0.7
   # pylab.figure()
   # pylab.contour(h,25)
   # pylab.axis('scaled')
   # pylab.xlim([0, h.shape[0]-1])
   # pylab.ylim([0, h.shape[0]-1])
   # pylab.show()
   # h[h < 0] = 0.0

   n = h.shape[0]
   dh = (numpy.arange(n) / float(n-1) - 0.5) / 2.0
   dh = dh ** 2 * n * SCALE
   for i in xrange(n):
      h[:,i] -= dh
      h[i,:] -= dh

   print 'Finished.'

   # build OpenGL list
   gllist = glGenLists(1)
   glNewList(gllist, GL_COMPILE)
   for i in xrange(1, n):
      glBegin(GL_QUAD_STRIP)
      for j in xrange(n):
	 glColor3f(r[j,i-1], g[j,i-1], b[j,i-1])
         glVertex3f((i-1) * SCALE, j * SCALE, h[j,i-1])
	 glColor3f(r[j,i], g[j,i], b[j,i])
         glVertex3f(i * SCALE, j * SCALE, h[j,i])
      glEnd()
   glEndList()

   glViewport(0, 0, 1200, 800)
   glShadeModel(GL_SMOOTH)
   # glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
   glEnable(GL_DEPTH_TEST)

   # atmosphere fog
   glClearColor(0.22, 0.45, 0.9, 1.0)
   glFogi(GL_FOG_MODE, GL_EXP)
   glFogfv(GL_FOG_COLOR, [0.22, 0.45, 0.9, 1.0])
   glFogf(GL_FOG_DENSITY, 0.025)
   glEnable(GL_FOG)

   glMatrixMode(GL_PROJECTION)
   glLoadIdentity()
   gluPerspective(45.0, 1.5, 0.1, 1000.0)
   glMatrixMode(GL_MODELVIEW)

   return gllist, h, dh

def drawScene():
   global dh
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
   glLoadIdentity()
   # calculate view orientation for greater height
   sight = z_eye * 4
   if sight > 128: sight = 128
   # calculate head up direction
   x, y = x_eye - h.shape[0]/2, y_eye - h.shape[0]/2
   dh_local = sum(dh[[x_eye, y_eye]]) * 2
   x_up, y_up = dh_local * x / float(x**2 + y**2) / SCALE, \
		dh_local * y / float(x**2 + y**2) / SCALE
   h_diff = (x_up * dx + y_up * dy) * sight
   gluLookAt(x_eye * SCALE, y_eye * SCALE, h[y_eye, x_eye] + z_eye, \
	     x_eye * SCALE + dx*sight, y_eye * SCALE + dy*sight, \
	     h[y_eye, x_eye] - h_diff, x_up, y_up, 1.0)
   glScalef(1.0, 1.0, 1.0)
   glCallList(gllist)
   glutSwapBuffers()

def keyPressed(key, x, y):
   global x_eye, y_eye, z_eye, dx, dy, fog
   if key == 'x':
      glutDestroyWindow(window)
      sys.exit()
   elif key == 'w':
      x_eye += dx
      y_eye += dy
   elif key == 's':
      x_eye -= dx
      y_eye -= dy
   elif key == 'a':
      x_eye -= dy
      y_eye += dx
   elif key == 'd':
      x_eye += dy
      y_eye -= dx
   elif key == 'j':
      z_eye += 0.25
   elif key == 'k':
      z_eye -= 0.25
      if z_eye < 0.25: z_eye = 0.25
   elif key == 'h':
      tmp = dy
      dy = dx
      dx = -tmp
   elif key == 'l':
      tmp = dy
      dy = -dx
      dx = tmp
   elif key == 'f':
      if fog:
	 glDisable(GL_FOG)
	 fog = False
      else:
	 glEnable(GL_FOG)
	 fog = True
   drawScene()

x_eye, y_eye, z_eye = 256, 128, 0.25
dx, dy = 0, 1
fog = True

glutInit(())
glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
glutInitWindowSize(1200, 800)
window = glutCreateWindow('continent')

glutDisplayFunc(drawScene)
glutKeyboardFunc(keyPressed)
gllist, h, dh = initScene()
glutMainLoop()

