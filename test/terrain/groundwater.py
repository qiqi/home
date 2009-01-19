import pylab
import numpy
from numpy import logical_or, logical_and, logical_not, zeros, random, \
                  maximum, minimum, unravel_index

from continent import continent

# def flowpath(h_node):
h_node = continent(4,5,0.465,8)    # penesulas

# calculate in which way water flows on heightmap
h = (h_node[1:,1:] + h_node[:-1,1:] + h_node[1:,:-1] + h_node[:-1,:-1]) / 4.0
sea_node = (h_node == 0)
sea = logical_or( logical_or(sea_node[1:,1:], sea_node[:-1,1:]), \
                  logical_or(sea_node[1:,:-1], sea_node[:-1,:-1]) )
# for i in range(1):
#    h = (h[1::2,1::2] + h[0::2,1::2] + h[1::2,0::2] + h[0::2,0::2]) / 4.0
#    sea = logical_or(logical_or(sea[1::2,1::2], sea[0::2,1::2]), \
#                     logical_or(sea[1::2,0::2], sea[0::2,0::2]))
h[sea] = 0.0

# dir = 0=+x, 1=-x, 2=+y, 3=-y, -1=sink (sea), -2=undetermined
direc = zeros(h.shape) - 2
direc[sea] = -1

# The following section calculates flow direction.
# Looping through nodes / cells are avoided by vectorization, for the
# performance of python loops is inferior.
assert h.shape[0] == h.shape[1]
n = h.shape[0]
indexi = [range(2,n), range(0,n-2), range(1,n-1), range(1,n-1)]
indexj = [range(1,n-1), range(1,n-1), range(2,n), range(0,n-2)]
# find out which cells are higher than surounding
h_ind = zeros([n, n, 4], dtype=bool)       # height indicator
for i in range(4):
   h_ind[1:n-1,1:n-1,i] = (h[1:n-1, 1:n-1] >= h[indexi[i],:][:,indexj[i]])

# cells that have only one descent direction
single = logical_and(h_ind.sum(2) == 1, direc == -2)
for i in range(4):
   direc[logical_and(single, h_ind[:,:,i])] = i

# randomize flow direction for the cells with multiple descent directions
multi = logical_and(h_ind.sum(2) > 1, direc == -2)
# based on height difference
h_diff = zeros([n, n, 4], dtype=float)
for i in range(4):
   h_diff[1:n-1,1:n-1,i] = h[1:n-1, 1:n-1] - h[indexi[i],:][:,indexj[i]]
h_diff[h_diff < 0.0] = 0.0
# normalize to sum to 1
h_diff[h_ind] += 0.0001
diff_sum = h_diff.sum(2)
for i in range(4):
   h_diff[:,:,i][multi] /= diff_sum[multi]
for i in range(1,4):
   h_diff[:,:,i] += h_diff[:,:,i-1]
# random select
select = random.random([n, n])
s_ind = zeros([n, n], dtype=int)       # random selection indicator
for i in range(3):
   s_ind[select > h_diff[:,:,i]] = i + 1
direc[multi] = s_ind[multi]

# group cells based on final destination
dest_i, dest_j = zeros([n, n], dtype=int), zeros([n, n], dtype=int)
for i in range(n):
   dest_i[i,:], dest_j[:,i] = i, i
i0, j0 = dest_i.copy(), dest_j.copy()
dir_next = direc[dest_i, dest_j]
while (dir_next >= 0).sum() > 0:
   print 'a'
   dest_i[dir_next == 0] += 1
   dest_i[dir_next == 1] -= 1
   dest_j[dir_next == 2] += 1
   dest_j[dir_next == 3] -= 1
   i_jmp = dest_i[dest_i, dest_j]
   j_jmp = dest_j[dest_i, dest_j]
   dest_i = i_jmp
   dest_j = j_jmp
   dir_next = direc[dest_i, dest_j]
dest_h = h[dest_i, dest_j]         # height of final destination
dest_i[dest_h == 0.0] = 0
dest_j[dest_h == 0.0] = 0

igorge = zeros([n-2, n-2], dtype = int)
hgorge = zeros([n-2, n-2], dtype = float) + 100.0
# for each point choose best direction to form gorge
for i in range(4):
   dest_h_nbr = dest_h[indexi[i],:][:,indexj[i]]
   idir = (dest_h[1:n-1,1:n-1] > dest_h_nbr)
   hnew = maximum(h[1:n-1,1:n-1], h[indexi[i],:][:,indexj[i]])
   hnew[logical_not(idir)] = 100.0
   igorge[hnew < hgorge] = i
   hgorge = minimum(hgorge, hnew)
# make them full size
tmp = (igorge, hgorge)
igorge = zeros([n, n], dtype = int)
hgorge = zeros([n, n], dtype = float) + 100.0
igorge[1:n-1,1:n-1], hgorge[1:n-1,1:n-1] = tmp

dead = logical_and(h_ind.sum(2) == 0, direc == -2)
dead_list, gorge_list = [], []
h_tmp = zeros([n, n])
ind = 0
for inj in set((dest_i * n + dest_j).flat):
   i, j = inj / n, inj % n
   print ind
   ind += 1
   if not i == j == 0:
      assert dead[i,j]
      dead_list.append((i, j))
      h_tmp[:,:] = 100.0
      h_tmp[dest_h == h[i,j]] = hgorge[dest_h == h[i,j]]
      ig, jg = unravel_index(h_tmp.argmin(), [n,n])







