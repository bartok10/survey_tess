# coding: utf-8
import numpy as np
randoms=np.loadtxt('/home/jarmijo/data/SDSS/random_SDSSDR7.dat')
rz = randoms[:,2]
rra = randoms[:,0]
rdec = randoms[:,1]
# %load 8
NSIDE=512
rrdec=np.radians(rdec)+np.pi/2.;rrra=np.radians(rra);
ramin=min(rrra)
ramax=max(rrra)
decmin = min(rrdec)
decmax = max(rrdec)
ranmask = np.loadtxt('data/SDSS/ranmaskDR7.dat')
zmin=0.0;zmax=0.05
edgemask=[]
ranRADecZ=[]
for i in range(len(randoms)):
    rara = np.random.uniform(ramin,ramax)
    radec = np.random.uniform(decmin,decmax)
    raz = np.random.uniform(zmin,zmax)
    p=hp.ang2pix(NSIDE,theta=radec,phi=rara,nest=True)
    if ranmask[p] == 0.0: 
        edgemask.append(p)
        ranRADecZ.append([rara,radec,raz])     
    
for i in range(len(randoms)/10):
    rara = np.random.uniform(ramin,ramax)
    radec = np.random.uniform(decmin,decmax)
    raz = np.random.uniform(0.0495,0.0505)
    p=hp.ang2pix(NSIDE,theta=radec,phi=rara,nest=True)
    if ranmask[p] == 1.0: tapa.append([rara,radec,raz])
        
