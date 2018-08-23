# coding: utf-8
# %load python/gen_mask2.py
import numpy as np
import healpy as hp
import sys,os
#dir1 = "/home/jarmijo/test_st/"
dir1 = "/home/joaquin/test_st/"
#print "Read random sample with the survey geometry... \n"
#ran_cat = '/home/jarmijo/SDSS/mock/random_SDSSDR7.dat' #sys.argv[1]
#ran_cat = '/home/joaquin/Documents/SDSS/mock/random_SDSSDR7.dat'
#randoms=np.loadtxt(ran_cat)
#rz = randoms[:,2]
#rra = randoms[:,0]
#rdec = randoms[:,1]
#
NSIDE=512
#rrdec=np.radians(rdec)+np.pi/2.;rrra=np.radians(rra); dec = 0 .. pi
#ramin=min(rrra);ramax=max(rrra);decmin = min(rrdec);decmax = max(rrdec) # limits from the mask
ramin = 1.916495; ramax = 4.566342; decmin = 1.505382; decmax = 2.797541;
rmsk = '/home/joaquin/Documents/SDSS/mock/ranmaskDR7.dat'
rana = np.loadtxt(rmsk) #mask with pixels filled by the randoms galaxies
ranmask = np.zeros(12*NSIDE**2)
ranmask[rana.astype(int)] = 1.
### select limits in redshift ###
print "select redshift limits... \n"
zmin=0.0;zmax=0.2 #zmin = sys.argv[2]; zmax= sys.argv[3]
nd = 0.558323 #mean number of galaxies per pixel in the total sample
edgemask=[]
ranRADecZ=[]
nl = 0.01 #new limit for edge mask
#
tapa1=[]
zl = 0.001 # thickness of each cap... it must be very thin

for i in range(int(5e5)):  # (N = np*NP) number density of randoms (10*np) comes here
    rara = np.random.uniform(ramin - nl, ramax + nl)
    radec = np.random.uniform(decmin - nl, decmax + nl)  # RaDec for random
    raz = np.random.uniform(zmin, zmax)  # uniform redshift distribution
    p = hp.ang2pix(NSIDE, theta=radec, phi=rara, nest=True)  # ang to pixel for random
    if ranmask[p] == 0.0:  ### if new random is outside the previous mask is an edge mock galaxy
        edgemask.append(p)
        ranRADecZ.append(
            [rara, radec, raz])  #save edge mock galaxynP = list(set(edgemask)) #number of pixels in the edge-mask
#
for i in range(int(5e4)): #different N (same number density)
    rara = np.random.uniform(ramin,ramax)
    radec = np.random.uniform(decmin,decmax)
    raz1 = np.random.uniform(zmax-zl,zmax+zl) # zmin needs to be higher than the thickne
    p=hp.ang2pix(NSIDE,theta=radec,phi=rara,nest=True)
    if ranmask[p] == 1.0: tapa1.append([rara,radec,raz1])

nP = list(set(edgemask)) #number of pixels in the edge-mask
ndp = float(len(ranRADecZ))/float(len(nP))
den_w = ndp / nd
print "wrapping is " + "%.3f" % den_w + " times more dense than the sample."
# save edges and caps (wrapping) RaDec, z
#np.savetxt(dir1+'mock_cap_b.txt',np.array(tapa1))
np.savetxt(dir1+'mock_cap_u.txt',np.array(tapa1))
np.savetxt(dir1+'mock_edges.txt',np.array(ranRADecZ))
# save edgmask pixels
np.savetxt(dir1+'edgemaskDR7.dat',np.array(edgemask))
