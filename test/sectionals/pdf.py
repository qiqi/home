import os
import math

import Image
from reportlab.pdfgen import canvas
from reportlab.lib import pagesizes

versions = ['42']

for version in versions:
  for chart in os.listdir(version):
    chartdir = os.path.join(version, chart)
    assert os.path.isdir(chartdir)
    print '    making pdf for chart %s...' % chartdir

    MARGIN = 25.0
    OVERLAP = 0.1
    
    im = Image.open( os.path.join(chartdir, 'chart.jpg') )
    page_size = pagesizes.letter[::-1]
    ratio = 2.60740740741
    content_size = tuple( x - 2*MARGIN for x in page_size )
    image_size = tuple( int(round(x * ratio)) for x in content_size )
    x0, y0 = 1000, 600
    nx = int(math.ceil((im.size[0] - x0) / float(image_size[0])))
    ny = int(math.ceil((im.size[1] - y0) / float(image_size[1])))
    dx = (im.size[0] - x0 - image_size[0]) / (nx - 1)
    dy = (im.size[1] - y0 - image_size[1]) / (ny - 1)
    
    c = canvas.Canvas(os.path.join(chartdir, 'chart.pdf'), pagesize=page_size)
    for ix in range(nx):
      for iy in range(ny):
        x = x0 + dx * ix
        y = y0 + dy * iy
        crop = im.crop((x, y, x + image_size[0], y + image_size[1]))
        crop_filename = 'crop%f%f.jpg' % (x, y)
        crop.save(crop_filename)
        c.drawImage(crop_filename, MARGIN, MARGIN, *content_size)
        os.remove(crop_filename)
        c.showPage()
    c.save()
