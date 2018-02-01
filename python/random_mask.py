# coding: utf-8
import healpy as hp
import numpy as np
pixels = np.loadtxt('/home/jarmijo/data/SDSS/maskhealpix.dr7.spec')
pixels -= 1 # to start from 0
NSIDE = 512
dim = 12*NSIDE**2 #dim mask
mask = np.zeros(dim)
mask[pixels.astype(int)] = 1
Npix = len(np.where(mask==1.0)[0])
pix1 = np.where(mask==1.0)[0] #pixeles con galaxias
NSIDE=512
pix0 = [] #pixeles sin galaxias 0 que tengan vecinos 1s
for i in pix1:
    vec = hp.get_all_neighbours(nside=512,theta=i,phi=None,nest=True)
    for j in vec:
        if mask[j]==0.0: pix0.append(j)
       
maskr = np.zeros(dim)#donde deberian ir randoms
pix0 = np.array(pix0)
maskr[pix0.astype(int)] = 1
#hp.mollview(map=maskr,nest=True,rot=(0.,180.),min=0,max=1,return_projected_map=True,title='DR7 masked',)
