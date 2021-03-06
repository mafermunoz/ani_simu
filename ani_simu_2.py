

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
dipole_dist=np.load('dipole_dist_events_ra_dec_ipix.npy')
##Sky with dipole_dist*10000*nevents/(Npix)
#Npix=healpy.nside2npix(NSIDE)
#dipole_dist=dipole_dist*100*(nevents/Npix) ## This is a healpy map
#sky=sky.astype(int)
##Create the list of tracks
#list_cr_ra=np.empty(sum(sky))
#list_cr_dec=np.empty(sum(sky))
#list_cr_npix=np.empty(sum(sky))
#list_cr_ra=np.array([])
#list_cr_dec=np.array([])
#list_cr_npix=np.array([])
#cpointing=0
ctracks=0

#print (np.sum(sky))

#print('here')
#for i in range (Npix):
#    ra,dec=healpy.pix2ang(nside=NSIDE,ipix=i)
#    dummy=np.ones(sky[i])
#    list_cr_ra=np.append(list_cr_ra,dummy*ra)
#    list_cr_dec=np.append(list_cr_dec,dummy*dec)
#    list_cr_npix=np.append(list_cr_npix,dummy*i)


#list=np.stack((list_cr_ra,list_cr_dec,list_cr_npix),axis=1)
len_list=len(dipole_dist)

print(len_list)
##The FOV of our simulated isntrument
fov=60

theta=np.empty(nevents)
phi=np.empty(nevents)
ra=np.empty(nevents)
dec=np.empty(nevents)

for i,x in enumerate(poisson_dist):
    if (i%10000)==0:
        print(i)
    if(i>0 and i%10000==0):
        np.savez('simu_info_ani3'+str(i),ra=ra,dec=dec,theta=theta,phi=phi)
    if(x==0):continue

    a=orbit[i]
    v=healpy.ang2vec(np.rad2deg(a[0]),np.rad2deg(a[1]),lonlat=True)
    search_radius=healpy.query_disc(NSIDE,v,np.deg2rad(fov))
    #print(search_radius.shape)

    c2=SkyCoord(np.rad2deg(a[0])*u.degree,np.rad2deg(a[1])*u.degree,frame='fk5')
    #print(c2)
    #print(search_radius)
    for j in range (x): ## Loops over the given amount of measurements per second
        dummy2=np.random.randint(len_list)
        aa=search_radius[search_radius==dipole_dist[dummy2,2]]
        #print(aa)
        while (aa.size>0)==False:
            del aa
            dummy2=np.random.randint(len_list)
            #print('stuck')
            #print(dipole_dist[dummy2,2])
            aa=search_radius[search_radius==dipole_dist[dummy2,2]]

        c1=SkyCoord(dipole_dist[dummy2,0]*u.degree,dipole_dist[dummy2,1]*u.degree,frame='fk5')
        #print(c1)
        offset=c2.spherical_offsets_to(c1)
        #print(offset)
        #print(np.where(search_radius==dipole_dist[dummy2,2]))
        ra[ctracks]=dipole_dist[dummy2,0]
        dec[ctracks]=dipole_dist[dummy2,1]
        theta[ctracks]=offset[0].degree
        phi[ctracks]=offset[1].degree
        dipole_dist[dummy2,2]=-55
        #print(ctracks)
        ctracks=ctracks+1
        del c1,offset

    del search_radius

np.savez('simu_info_ani3_all.npz',ra=ra,dec=dec,theta=theta,phi=phi)
