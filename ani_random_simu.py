import numpy as np
import healpy as hp
import sys

def thetaphi2radec(theta, phi, ra, dec, rotator=None):
    '''
    Transforms the satelite-based coordinates theta phi
    to the earth-based coordinates crra crdec
    (ra and dec of the cosmic ray)
    using the satelite pointing ra dec
    '''
    rot = [0., 90. - dec, -ra - 90.]

    if rotator is None:
        r = hp.Rotator(rot=rot, eulertype='X')
    else:
        r = rotator
        r._rots = [np.radians(rot)]
        r._update_matrix()

    theta = 90. - theta
#     crra, crdec = r(theta, phi, lonlat=True)
    crra, crdec = r(phi, theta, lonlat=True)

    if rotator is None:
        return crra, crdec, r
    return crra, crdec

def main(nmap,njob):
    nmap=int(nmap)##From 1 to 10
    njob=int(njob)
    poisson_dist=np.load('possion_split_all.npy')#From  file with the random tracks according to the poisson distribution
    start=0
    start_track=0
    for i in range(njob-1):
        start=start+len(poisson_dist[nmap-1][i])
        start_track=start_track+poisson_dist[nmap-1][i].sum()
    end=start+len(poisson_dist[nmap-1][njob-1])
    print(start)

    poisson_dist=poisson_dist[nmap-1][njob-1]
    track=np.load('tracks_shuffled_'+str(nmap-1)+'.npy')
    sat_pointing=np.load('DAMPE_2A_OBS_2016averge_pos_per_second.npy')##From Second based information
    nTracks=len(track)
    ctrack=start_track
    ra=np.empty(np.sum(poisson_dist))
    dec=np.empty(np.sum(poisson_dist))
    #dec=np.empty(nTracks)
    #ra=np.empty(nTracks)
    #dec=np.empty(nTracks)

    for i,x in enumerate (poisson_dist): ##Loop over the elements of the poison distribution
        if(i%100000==0):
            print(i)
        if x==0:
            continue
        for j in range (x):## For each of the tracks that the possion distribution says to us

            if(ctrack<nTracks): ## Counter to check we dont repeat tracks
                a=sat_pointing[start+i]
                ra[ctrack],dec[ctrack],r=thetaphi2radec(np.rad2deg(track[ctrack,1]),np.rad2deg(track[ctrack,2]),np.rad2deg(a[0]),np.rad2deg(a[1]))
                ctrack=ctrack+1
            else:
                #option to use the full statistics comming from the poisson Distribution
                a=sat_pointing[start+i]
                dd=np.random.randint(nTracks, size=1)
                ra[ctrack],dec[ctrack],r=thetaphi2radec(np.rad2deg(track[dd[0],1]),np.rad2deg(track[dd[0],2]),np.rad2deg(a[0]),np.rad2deg(a[1]))
                ctrack=ctrack+1
                break


    #ra=np.array(ra)
    #dec=np.array(dec)
    np.savez("../Map_3_"+str(nmap-1)+"_"+str(njob-1),ra=ra,dec=dec)




if __name__ == '__main__':
        main(sys.argv[1],sys.argv[2])
