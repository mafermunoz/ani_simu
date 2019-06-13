

import healpy
import numpy as np
import math
import operator
import sys
import glob
import astropy
from astropy.coordinates import SkyCoord  # High-level coordinates
import astropy.units as u


##Real Orbit for DAMPE for one year
orbit=np.load("/beegfs/dampe/users/mmunozsa/ani_random_average/ani_avg_method/DAMPE_2A_OBS_2016averge_pos_per_second.npy")
NSIDE=16
##Poissonian distrbution with the real data information
poisson_dist=np.load('../possion_mc_11.npy')
#nevents=sum of values from the possonian distribution.
nevents=np.sum(poisson_dist)
#npos_sky=np.where(poisson_dist>0)
##Use/define a dipole  of strenght 0.1 +1
dipole_dist=np.load('../dipole_distribution_nbins16_strenght_1.npy')
##Sky with dipole_dist*10000*nevents/(Npix)
Npix=healpy.nside2npix(NSIDE)
sky=dipole_dist*10000*(nevents/Npix) ## This is a healpy map
sky=sky.astype(int)
##Create the list of tracks
#list_cr_ra=np.empty(sum(sky))
#list_cr_dec=np.empty(sum(sky))
#list_cr_npix=np.empty(sum(sky))
list_cr_ra=np.array([])
list_cr_dec=np.array([])
list_cr_npix=np.array([])
cpointing=0
ctracks=0

print (np.sum(sky))

for i in range (Npix):
    ra,dec=healpy.pix2ang(nside=NSIDE,ipix=i)
    dummy=np.ones(sky[i])
    list_cr_ra=np.append(list_cr_ra,dummy*ra)
    list_cr_dec=np.append(list_cr_dec,dummy*dec)
    list_cr_npix=np.append(list_cr_npix,dummy*i)


list=np.stack((list_cr_ra,list_cr_dec,list_cr_npix),axis=1)
len_list=len(list)
print('len_list')
##The FOV of our simulated isntrument
fov=60

theta=np.empty(nevents)
phi=np.empty(nevents)
ra=np.empty(nevents)
dec=np.empty(nevents)


for i,x in enumerate(poisson_dist[:100]):
    ##x tells us the number of tracks at  a given dampe position
    if (x==0): continue

    a=orbit[i]
    v=healpy.ang2vec(a[0],a[1])
    search_radius=healpy.query_disc(NSIDE,v,np.deg2rad(fov),inclusive=True)
    c2=SkyCoord(a[0]*u.rad,a[1]*u.rad,frame='fk5')


    for j in range (x): ## Loops over the given amount of measurements per second
        dummy2=np.random.randint(len_list)

        while (list[dummy2,2] in search_radius)==False:
            dummy2=np.random.randint(len_list)

        c1=SkyCoord(list[dummy2,0]*u.rad,list[dummy2,1]*u.rad,frame='fk5')
        offset=c2.spherical_offsets_to(c1)
        ra[ctracks]=list[dummy2,0]
        dec[ctracks]=list[dummy2,1]
        theta[ctracks]=offset[0].rad
        phi[ctracks]=offset[1].rad
        list.pop(dummy2)
        len_list=len_list-1
        ctracks+=ctracks

np.savez('simu_info_ani.npz',ra=ra,dec=dec,theta=theta,phi=phi)