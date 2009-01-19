import os
import urllib

versions = ['42']
charts = map(lambda a: a.strip(), file('list.txt').readlines())

for version in versions:
  if os.path.exists(version):
    assert os.path.isdir(version)
  else:
    os.mkdir(version)
  for chart in charts:
    print 'getting chart ' + chart + '...'
    chartdir = os.path.join(version, chart)
    if os.path.exists(chartdir):
      assert os.path.isdir(chartdir)
    else:
      os.mkdir(chartdir)
    ni, nj = 0, 0
    for i in range(100):
      for j in range(100):
         url = 'http://skyvector.com/newtiles/%s/%s/tile_1_%d_%d.jpg' % \
               (version, chart, i, j)
         filename = os.path.join(chartdir, 'tile_1_%d_%d.jpg' % (i, j))
         urllib.urlretrieve(url, filename)
         if file(filename).readline().startswith('<!DOCTYPE HTML'):
           os.remove(filename)
           if j > 0:
             assert nj == 0 or nj == j
             nj = j
           break
      if j == 0:
        ni = i
        break
      print '    column %d obtained, totally %d tiles.' % (i, j)
    sizefile = os.path.join(chartdir, 'size.txt')
    file(sizefile, 'w').write('%d, %d' % (ni, nj))

