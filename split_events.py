import numpy as np
import math
import operator
import sys
import glob

orbit=np.load("/beegfs/dampe/users/mmunozsa/ani_random_average/ani_avg_method/DAMPE_2A_OBS_2016averge_pos_per_second.npy")
np.random.shuffle(orbit)
o_100=np.array_split(orbit,100)
print(o_100.shape())
np.save('../orbit_2016_j100',o_100)


poisson_dist=np.load('../possion_mc_11.npy')
p_100=np.array_split(poisson_dist,100)
print(p_100.shape())
np.save('../poison_2016_j100',p_100)

dipole_dist=np.load('dipole_dist_events_ra_dec_ipix.npy')
np.random.shuffle(dipole_dist)
d_100=np.array_split(dipole_dist)
print(dipole_dist.shape())
np.save('../dipole_dist_2016_j100',d_100)
