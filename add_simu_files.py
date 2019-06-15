import numpy as np
import math
import operator
import sys
import glob


file_path='data/*'
txt=glob.glob(file_path)
print(txt)
for i,file  in enumerate (txt):
    if i==0:
        f=np.load(file)
        ra=f['ra']
        dec=f['dec']
        theta=f['theta']
        phi=f['phi']
    else:
        dummy=np.load(file)
        ra=np.hstack([ra,dummy['ra']])
        dec=np.hstack([dec,dummy['dec']])
        theta=np.hstack([theta,dummy['theta']])
        phi=np.hstack([phi,dummy['phi']])

np.savez('ani_simu_all_all.npz',ra=ra,dec=dec,theta=theta,phi=phi)
