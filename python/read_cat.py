import numpy as np
import healpy as hp
dir1= '/home/jarmijo/test_st/'
cat = '/home/jarmijo/SDSS/data/sdssclean.dat'  #catalog to read
nmask = '/home/jarmijo/SDSS/mock/maskhealpix.dr7.spec' # SDSS mask to select the galaxies
zmin = 0.05; zmax = 0.1; Mrl = -20.4 ## dim2 as example
pos = np.loadtxt(cat,usecols=(0,1,2)) #Ra,Dec,z
Mr = np.loadtxt('/home/jarmijo/SDSS/data/sdss_corr.dat',usecols=(3,)) # Magnitude M_r to cut Main sample
pix = np.loadtxt(nmask)
NSIDE = 512; dim = 12*NSIDE**2;
mask = np.zeros(dim); pix -=1; 
mask[pix.astype(int)] = 1 #mask in array format

#reading each column (angle to radians)
z = pos[:,2]
ra = pos[:,0]
dec = pos[:,1]
ra_rad =np.radians(ra); dec_rad = np.radians(dec);dec_rad += np.pi/2. # from 0 to pi
#convert to healpix pixels
a = hp.pixelfunc.ang2pix(nside=NSIDE,theta=dec_rad,phi=ra_rad,nest=True)
gal_pix = np.sort(a); mgal = mask[gal_pix]
ra_m = ra_rad[mgal.astype(bool)];dec_m = dec_rad[mgal.astype(bool)]
z_m = z[mgal.astype(bool)]; Mr_m = Mr[mgal.astype(bool)]
# select subsample
ra_cat = ra_m[(z_m<zmax)&(z_m>zmin)&(Mr_m>Mrl)]; dec_cat = dec_m[(z_m<zmax)&(z_m>zmin)&(Mr_m>Mrl)];z_cat = z_m[(z_m<zmax)&(z_m>zmin)&(Mr_m>Mrl)]
#big patch
ra_b = ra_cat[(ra_cat>-1.3 + np.pi)&(ra_cat<1.55+np.pi)]; dec_b = dec_cat[(ra_cat>-1.3 + np.pi)&(ra_cat<1.55+np.pi)];z_b = z_cat[(ra_cat>-1.3 + np.pi)&(ra_cat<1.55+np.pi)]
ra_up = ra_b[(ra_b<1.2+np.pi)&(dec_b > 0.8 + np.pi/2.)]
dec_up = dec_b[(ra_b<1.2+np.pi)&(dec_b > 0.8 + np.pi/2.)]
z_up = z_b[(ra_b<1.2+np.pi)&(dec_b > 0.8 + np.pi/2.)]
#
ra_bot = ra_b[(ra_b<1.5+np.pi)&(dec_b < 0.8 + np.pi/2.)]
dec_bot = dec_b[(ra_b<1.5+np.pi)&(dec_b < 0.8 + np.pi/2.)]
z_bot = z_b[(ra_b<1.5+np.pi)&(dec_b < 0.8 + np.pi/2.)]
ra_p = np.concatenate([ra_up,ra_bot]);dec_p = np.concatenate([dec_up,dec_bot]);z_p = np.concatenate([z_up,z_bot])
S = np.array([ra_p,dec_p,z_p])
np.savetxt(dir1+'subsample.txt',S.T)
print "file saved in:" + dir1 + 'subsample.txt' 
print "ng:" + str(len(z_p)) + " np:" + str(len(pix)) + "\n"
dp = float(len(z_p))/float(len(pix))
print "mean number of galaxies per pixel in the sample:" + str(dp)
#ng = len(z_p); NP = len(pix);



