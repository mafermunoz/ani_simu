

import healpy
import numpy as np
import math
import operator
import sys
import glob

orbit=np.load("Orbit_dampe_2016.npy")
fake_events=np.load("fake_events.npy")


sep_max=np.deg2rad(60)
track=np.zeros((len(fake_events),2))
pos=np.zeros((len(fake_events),2))
for i,x in enumerate(fakes):
    if(i%1000==0):
        print(i)
    for j,y in enumerate(orbs):
        sep=healpy.rotator.angdist(x,y)
        if(sep<=sep_max):
                track[i]=abs(x-y)
                pos[i]=y
                orbs=np.delete(orbs,j,0)
                break
