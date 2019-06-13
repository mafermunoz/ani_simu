
import healpy
import numpy as np
import math
import operator
import sys
import glob
import astropy
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u

NSIDE=16
dipole_dist=np.load('../dipole_distribution_nbins16_strenght_1.npy')
##Sky with dipole_dist*10000*nevents/(Npix)
Npix=healpy.nside2npix(NSIDE)
nevents=5014862
sky=dipole_dist*100*(nevents/Npix) ## This is a healpy map
sky=sky.astype(int)
fake_events=np.zeros((sky.sum(),3))
count=0
for i in range (NPIX):
    ang=healpy.pix2ang(NSIDE,i,lonlat=True)
    print(ang)
    for j in range (sky[i]):
        a=count+j
        fake_events[a,:2]=ang
        fake_events[a,2]=i
    print(fake_events[count])
    count=count+simu_int[i]
np.save('dipole_dist_events_ra_dec_ipix.npy',fake_events)
