import os

import Image

versions = ['42']

for version in versions:
  for chart in os.listdir(version):
    chartdir = os.path.join(version, chart)
    assert os.path.isdir(chartdir)
    print '    concatenating chart %s...' % chartdir

    ni, nj = eval( file(os.path.join(chartdir, 'size.txt')).read() )
    
    # get total width and height
    width, height = [0], [0]
    for i in range(ni):
      filename = os.path.join(chartdir, 'tile_1_%d_%d.jpg' % (i, 0))
      width.append(width[-1] + Image.open(filename).size[1])
    for j in range(nj):
      filename = os.path.join(chartdir, 'tile_1_%d_%d.jpg' % (0, j))
      height.append(height[-1] + Image.open(filename).size[0])
    
    # concatenate image
    image = Image.new('RGB', (width[-1], height[-1]))
    for i in range(ni):
      for j in range(nj):
        filename = os.path.join(chartdir, 'tile_1_%d_%d.jpg' % (i, j))
        subimage = Image.open(filename)
        box = (width[i], height[j], width[i+1], height[j+1])
        image.paste(subimage, box)
    
    print '    saving concatenated image...'
    image.save(os.path.join(chartdir, 'chart.jpg'))
  
