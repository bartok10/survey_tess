# coding: utf-8
# %load python/RaDec2cartesian.py
import numpy as np
from astropy.cosmology import WMAP7 as cosmo
#dir1= '/home/jarmijo/test_st/'
dir1 = '/home/joaquin/test_st/'
gal_pos = np.loadtxt(dir1+'subsample.txt')
edges = np.loadtxt(dir1+'mock_edges.txt')
cap1 = np.loadtxt(dir1+'mock_cap_u.txt')
#cap2 = np.loadtxt(dir1+'mock_cap_b.txt')
#
c=299792.458;H0=68
dc_z = cosmo.comoving_distance(gal_pos[:,2]).to_value()
ra = gal_pos[:,0]; dec = gal_pos[:,1]
# se puede hacer lo mismo con rutinas  healpy
############# galaxias ############
xp = dc_z * np.cos(ra) * np.sin(dec)
yp = dc_z * np.sin(ra) * np.sin(dec)
zp = dc_z * np.cos(dec)
########### edges (and holes) ############
xe=cosmo.comoving_distance(edges[:,2]).to_value() * np.cos(edges[:,0]) * np.sin(edges[:,1])
ye=cosmo.comoving_distance(edges[:,2]).to_value() * np.sin(edges[:,0]) * np.sin(edges[:,1])
ze=cosmo.comoving_distance(edges[:,2]).to_value() * np.cos(edges[:,1])
########### cap ################
xc1=cosmo.comoving_distance(cap1[:,2]).to_value() * np.cos(cap1[:,0]) * np.sin(cap1[:,1])
yc1=cosmo.comoving_distance(cap1[:,2]).to_value() * np.sin(cap1[:,0]) * np.sin(cap1[:,1])
zc1=cosmo.comoving_distance(cap1[:,2]).to_value() * np.cos(cap1[:,1])
#
np.savetxt(dir1+'gal_pos_3d_com.txt',np.array([xp,yp,zp]).T)
p_file = open(dir1+'cat_zobov.pos','a')
p_file.write(str(len(xp)+len(xe)+len(xc1)) + ' ' + str(len(xp)) + '\n')
for i in range(len(xp)):
    p_file.write('%.8f' %  xp[i] + ' ' + '%.8f' %  yp[i] + ' ' + '%.8f' %  zp[i] + '\n')
for i in range(len(xe)):
    p_file.write('%.8f' %  xe[i] + ' ' + '%.8f' %  ye[i] + ' ' + '%.8f' %  ze[i] + '\n')
for i in range(len(xc1)):
    p_file.write('%.8f' %  xc1[i] + ' ' + '%.8f' %  yc1[i] + ' ' + '%.8f' %  zc1[i] + '\n')
    
p_file.close()
