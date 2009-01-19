from numpy import log, arange, cos, sin, sqrt, zeros, pi, random, power, amin

from continent import randimage

def temperature(h, lat0, lat1):
   # calculate temperature based on heightmap and latitude
   # return summer and winter temperature
   lat = arange(h.shape[1]) / float(h.shape[1] - 1) * (lat1 - lat0) + lat0
   lat = abs(lat) / 90.0
   assert lat.max() >= 0.0 and lat.max() <= 1.0
   # sea level temperature for each latitude (centigrade)
   # City       latitude          winter          summer  (assuming at sea lvl)
   # Hongkong     22.5             +10              37
   # Shanghai     30.0              +3              34
   # Beijing      40.0              -2              32
   # Haerbin      45.0              -8              28
   # Mosco        55.0             -13              23
   s = cos(lat * pi/2) * 40.0
   w = cos(sqrt(lat) * pi) * 30.0 + 10.0
   # import pylab; pylab.plot(s); pylab.plot(w)
   summer, winter = zeros(h.shape), zeros(h.shape)
   for i in range(h.shape[0]):
      summer[:,i], winter[:,i] = s, w
   # adjust for altitude (-13 degrees per km)
   summer -= h * 7.0
   winter -= h * 6.0
   reversed = summer < winter
   summer[reversed] = winter[reversed]
   return summer, winter

def precipitation(h, lat0, lat1):
   # calculate precipitation based on heightmap and latitude
   assert h.shape[0] == h.shape[1]
   n = h.shape[0] - 1
   niter = 0
   while n > 2:
       assert n % 4 == 0
       n /= 2
       niter += 1
   # generate course grid height map and roughness map
   dh = abs(h[1:,1:] - h[1:,:-1]) + abs(h[:-1,1:] - h[:-1,:-1]) + \
        abs(h[1:,1:] - h[:-1,1:]) + abs(h[1:,:-1] - h[:-1,:-1])
   hc = amin([h[1:,1:], h[1:,:-1], h[:-1,1:], h[:-1,:-1]], axis=0)
   nc = 6
   dx = 1
   for iter in range(nc, niter):
      dx *= 2
      dh = (dh[0::2,0::2] + dh[1::2,0::2] + dh[0::2,1::2] + dh[1::2,1::2]) / 4.0
      hc = (hc[0::2,0::2] + hc[1::2,0::2] + hc[0::2,1::2] + hc[1::2,1::2]) / 4.0
   dh = (dh - dh.min()) / (dh.max() - dh.min())
   # generate random wind direction
   nwind = 3
   wind_str, wind_dir = random.random(nwind), random.random(nwind)
   windx = wind_str * cos(wind_dir * 2*pi) / dx
   windy = wind_str * sin(wind_dir * 2*pi) / dx
   wind_season = random.random(nwind)
   # propagate vapor with wind
   precip_total = zeros(dh.shape)
   for iwind in range(nwind):
      # set up constant source condition: vapor in sea
      summer, winter = temperature(hc, lat0, lat1)
      temp = wind_season[iwind] * (winter - summer) + summer
      source = power(10.0, temp / 40.0 - 1.0)  # 0 degree is 10% of 40 degree...
      source[hc > 0] = 0.0     # no source on land
      # initialize density of water contained in air
      dt = 1.0
      precip_rate = 0.001  # rate of water in air forms precipitation
      rate = precip_rate * (1 + 50.0 * sqrt(dh))
      # taking consideration of barrier effect of terrain
      density = source / precip_rate
      # time iteration
      for istep in range(1000):
         # source and sink
         density += source * dt
         precip = density * rate
         density -= precip * dt
         # convection
         if windx[iwind] > 0.0:
            density[1:,:] -= windx[iwind] * (density[1:,:] - density[0:-1,:])
         else:
            density[0:-1,:] -= windx[iwind] * (density[1:,:] - density[0:-1,:])
         if windy[iwind] > 0.0:
            density[:,1:] -= windy[iwind] * (density[:,1:] - density[:,0:-1])
         else:
            density[:,0:-1] -= windy[iwind] * (density[:,1:] - density[:,0:-1])
      precip_total += precip
   # interpolate back
   for iter in range(nc, niter):
      p0 = precip_total
      precip_total = zeros([p0.shape[0] * 2, p0.shape[1] * 2])
      precip_total[0::2,0::2] = p0
      precip_total[1::2,0::2] = p0
      precip_total[0::2,1::2] = p0
      precip_total[1::2,1::2] = p0
      precip_total[1:-1,1:-1] += (-4 * precip_total[1:-1,1:-1] + \
                precip_total[:-2,1:-1] + precip_total[2:,1:-1] + \
                precip_total[1:-1,:-2] + precip_total[1:-1,2:]) / 8.0
   # from cell centered to mesh centered
   p0 = precip_total / 4.0
   precip_total = zeros([p0.shape[0] + 1, p0.shape[1] + 1])
   precip_total[1:,1:] += p0
   precip_total[:-1,1:] += p0
   precip_total[1:,:-1] += p0
   precip_total[:-1,:-1] += p0
   precip_total[[0,-1],:] *= 2.0
   precip_total[:,[0,-1]] *= 2.0
   # normalize
   assert precip_total.min() > 0
   precip_total /= precip_total.max()
   return precip_total

if __name__ == '__main__':
   import pylab
   from continent import continent

   h = continent(13,13,0.0,8)    # good plains
   
   pylab.figure()
   s, w = temperature(h*8, 10, 60)
   pylab.contour(h,[0],colors='k')
   pylab.contour(s,25,colors='r')
   pylab.contour(w,25,colors='b')
   pylab.axis('scaled')
   pylab.xlim([0, h.shape[0]-1])
   pylab.ylim([0, h.shape[0]-1])
   pylab.title('temperature')
   
   p = precipitation(h, 10, 60)
   
   pylab.figure()
   pylab.contour(h,[0],colors='k')
   pylab.contour(p,[0.4,0.2,0.1,0.05,0.025,0.01])
   pylab.axis('scaled')
   pylab.xlim([0, h.shape[0]-1])
   pylab.ylim([0, h.shape[0]-1])
   pylab.title('precipitation')
