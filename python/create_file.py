# coding: utf-8
# %load python/RaDec2cartesian.py
import numpy as np
gal_pos = np.loadtxt('/home/jarmijo/SDSS/data/dim1/pos_RaDecZ.dat')
edges = np.loadtxt('/home/jarmijo/SDSS/data/dim1/mock_edges.txt')
cap = np.loadtxt('/home/jarmijo/SDSS/data/dim1/mock_cap.txt')
ra_d = gal_pos[:,0];dec_d=gal_pos[:,1];z_d=gal_pos[:,2]
ra_b = ra_d[(ra_d>-1.3 + np.pi)&(ra_d<1.55+np.pi)]; dec_b = dec_d[(ra_d>-1.3 + np.pi)&(ra_d<1.55+np.pi)];z_b = z_d[(ra_d>-1.3 + np.pi)&(ra_d<1.55+np.pi)]
ra_up = ra_b[(ra_b<1.2+np.pi)&(dec_b > 0.8 + np.pi/2.)]
dec_up = dec_b[(ra_b<1.2+np.pi)&(dec_b > 0.8 + np.pi/2.)]
z_up = z_b[(ra_b<1.2+np.pi)&(dec_b > 0.8 + np.pi/2.)]
#
ra_bot = ra_b[(ra_b<1.5+np.pi)&(dec_b < 0.8 + np.pi/2.)]
dec_bot = dec_b[(ra_b<1.5+np.pi)&(dec_b < 0.8 + np.pi/2.)]
z_bot = z_b[(ra_b<1.5+np.pi)&(dec_b < 0.8 + np.pi/2.)]
ra_p = np.concatenate([ra_up,ra_bot]);dec_p = np.concatenate([dec_up,dec_bot]);z_p = np.concatenate([z_up,z_bot])
c=299792.458;H0=68
# se puede hacer lo mismo con rutinas  healpy
############# galaxias ############
xp = c*z_p/H0 * np.cos(ra_p) * np.cos(dec_p)
yp = c*z_p/H0 * np.sin(ra_p) * np.cos(dec_p)
zp = c*z_p/H0 * np.sin(dec_p)
########### edges (and holes) ############
xe=c*edges[:,2]/H0 * np.cos(edges[:,0]) * np.cos(edges[:,1])
ye=c*edges[:,2]/H0 * np.sin(edges[:,0]) * np.cos(edges[:,1])
ze=c*edges[:,2]/H0 * np.sin(edges[:,1])
########### cap ################
xc=c*cap[:,2]/H0 * np.cos(cap[:,0]) * np.cos(cap[:,1])
yc=c*cap[:,2]/H0 * np.sin(cap[:,0]) * np.cos(cap[:,1])
zc=c*cap[:,2]/H0 * np.sin(cap[:,1])

p_file = open('/home/jarmijo/SDSS/data/dim1/test.pos','a')
p_file.write(str(len(xp)) + ' ' + str(len(xe) + len(xc)) + '\n')
for i in range(len(xp)):
    p_file.write('\t' + '%.5f' %  xp[i] + '\t' + '%.5f' %  yp[i] + '\t' + '%.5f' %  zp[i] + '\n')
for i in range(len(xe)):
    p_file.write('\t' + '%.5f' %  xe[i] + '\t' + '%.5f' %  ye[i] + '\t' + '%.5f' %  ze[i] + '\n')
for i in range(len(xc)):
    p_file.write('\t' + '%.5f' %  xc[i] + '\t' + '%.5f' %  yc[i] + '\t' + '%.5f' %  zc[i] + '\n')
    
p_file.close()
