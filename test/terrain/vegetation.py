from numpy import zeros, logical_and

def vegetation(h, dh, summer, winter, precip):
   # colormap
   # r = continent(2,4,None,8)
   # g = continent(3,8,None,8)
   # b = continent(0,7,None,8)
   r, g, b, t = zeros(h.shape), zeros(h.shape), zeros(h.shape), zeros(h.shape)
   
   # desert
   r[:] = 0.9 - (dh % 0.0001) * 2000.0
   g[:] = 0.9 - (dh % 0.0001) * 2000.0
   b[:] = 0.2 + (dh % 0.0001) * 3000.0
   t[:] = 100
   # tundra
   tundra = logical_and(precip > 0.01, summer > 5)
   r[tundra] = 0.5
   g[tundra] = 0.6 - (dh[tundra] % 0.01) * 10
   b[tundra] = 0.2
   t[tundra] = 110
   # cold grassland
   grass = logical_and(precip > 0.025, winter + summer > -5)
   r[grass] = 0.6
   g[grass] = 0.8 - dh[grass] * 0.2
   b[grass] = 0.3
   t[grass] = 120
   # temporal grassland
   grass2 = logical_and(precip > 0.025, winter + summer > 15)
   r[grass2] = 0.4
   g[grass2] = 0.8 - dh[grass2] * 0.1
   b[grass2] = 0.3
   t[grass2] = 220
   # conifer forest
   conifer = logical_and(precip > 0.08, summer > 5)
   conifer = logical_and(conifer, summer + winter > 0)
   r[conifer] = 0.0
   g[conifer] = 0.4
   b[conifer] = 0.1
   t[grass2] = 150
   # subtrop grassland
   grass3 = logical_and(precip > 0.025, winter + summer > 25)
   r[grass3] = 0.4
   g[grass3] = 0.6 - dh[grass3] * 0.1
   b[grass3] = 0.3
   t[grass2] = 320
   # broadleave forest
   broad = logical_and(precip > 0.12, summer > 30)
   r[broad] = 0.3
   g[broad] = 0.5
   b[broad] = 0.0
   t[grass2] = 350
   # trop desert
   desert = logical_and(precip < 0.2, winter + summer > 40)
   r[desert] = 0.9 - (dh[desert] % 0.0001) * 2000.0
   g[desert] = 0.9 - (dh[desert] % 0.0001) * 2000.0
   b[desert] = 0.2 + (dh[desert] % 0.0001) * 3000.0
   t[grass2] = 400
   # trop grassland
   grass4 = logical_and(precip > 0.2, winter > -5)
   grass4 = logical_and(grass4, dh < 0.01)
   r[grass4] = 0.5
   g[grass4] = 0.6 - dh[grass4] * 0.1
   b[grass4] = 0.1
   t[grass2] = 420
   # tropical rainforest
   trop_forest = logical_and(precip > 0.4, winter > 13)
   r[trop_forest] = 0.0
   g[trop_forest] = 0.5 - dh[trop_forest] * 0.1
   b[trop_forest] = 0.1
   t[grass2] = 450
   # snow cover
   snow = summer <= 0
   r[snow] = 1.0
   g[snow] = 1.0
   b[snow] = 1.0
   t[grass2] = 0
   # sea
   r[h == 0] = 0.0
   g[h == 0] = 0.1
   b[h == 0] = 0.6
   t[grass2] = 900

   return r, g, b, t
