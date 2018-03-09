# coding: utf-8
# %load python/RaDec2cartesian.py
import numpy as np
dir1= '/home/jarmijo/test_st/'
gal_pos = np.loadtxt(dir1+'subsample.txt')
edges = np.loadtxt(dir1+'mock_edges.txt')
cap1 = np.loadtxt(dir1+'mock_cap_u.txt')
cap2 = np.loadtxt(dir1+'mock_cap_b.txt')
#
c=299792.458;H0=68
# se puede hacer lo mismo con rutinas  healpy
############# galaxias ############
xp = c*gal_pos[:,2]/H0 * np.cos(gal_pos[:,0]) * np.cos(gal_pos[:,1])
yp = c*gal_pos[:,2]/H0 * np.sin(gal_pos[:,0]) * np.cos(gal_pos[:,1])
zp = c*gal_pos[:,2]/H0 * np.sin(gal_pos[:,1])
########### edges (and holes) ############
xe=c*edges[:,2]/H0 * np.cos(edges[:,0]) * np.cos(edges[:,1])
ye=c*edges[:,2]/H0 * np.sin(edges[:,0]) * np.cos(edges[:,1])
ze=c*edges[:,2]/H0 * np.sin(edges[:,1])
########### cap ################
xc1=c*cap1[:,2]/H0 * np.cos(cap1[:,0]) * np.cos(cap1[:,1])
yc1=c*cap1[:,2]/H0 * np.sin(cap1[:,0]) * np.cos(cap1[:,1])
zc1=c*cap1[:,2]/H0 * np.sin(cap1[:,1])
#
xc2=c*cap2[:,2]/H0 * np.cos(cap2[:,0]) * np.cos(cap2[:,1])
yc2=c*cap2[:,2]/H0 * np.sin(cap2[:,0]) * np.cos(cap2[:,1])
zc2=c*cap2[:,2]/H0 * np.sin(cap2[:,1])

#xp -= (xe.min()- 0.1);yp -= (ye.min() - 0.1)
#xc1 -= (xe.min()- 0.1);yc1 -= (ye.min() - 0.1)
#xc2 -= (xe.min()- 0.1);yc2 -= (ye.min() - 0.1)
#xe -= (xe.min()- 0.1);ye -= (ye.min() - 0.1)
np.savetxt(dir1+'gal_pos.txt',np.array([xp,yp,zp]).T)
p_file = open(dir1+'cat_zobov.pos','a')
p_file.write(str(len(xp)+len(xe)+len(xc1)+len(xc2)) + ' ' + str(len(xp)) + '\n')
for i in range(len(xp)):
    p_file.write('%.8f' %  xp[i] + ' ' + '%.8f' %  yp[i] + ' ' + '%.8f' %  zp[i] + '\n')
for i in range(len(xe)):
    p_file.write('%.8f' %  xe[i] + ' ' + '%.8f' %  ye[i] + ' ' + '%.8f' %  ze[i] + '\n')
for i in range(len(xc1)):
    p_file.write('%.8f' %  xc1[i] + ' ' + '%.8f' %  yc1[i] + ' ' + '%.8f' %  zc1[i] + '\n')
for i in range(len(xc2)):
    p_file.write('%.8f' %  xc2[i] + ' ' + '%.8f' %  yc2[i] + ' ' + '%.8f' %  zc2[i] + '\n')

    
p_file.close()
