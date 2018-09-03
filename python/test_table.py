# coding: utf-8
import numpy as np
from surveytess import  ZobovTess,voronoi_properties
from scipy.spatial import Voronoi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
#
#load files
# dir1 = '/home/jarmijo/CHUVIS/'
dir1 = '/home/ntejos/python/survey_tess/CHUVIS/'
gal = np.loadtxt(dir1+'gal_XYZ.dat')
gal_rdz = np.loadtxt(dir1+'subsample.dat')
vol = np.loadtxt(dir1+'cat_zobov.ascii.vol',skiprows=1)
zones = np.loadtxt(dir1+'cat_zobov.zone',skiprows=1)
f = open(dir1+'cat_zobov.void','r')
af  = f.read().splitlines()
f.close()
adj_filename = dir1+"cat_zobov.ascii.adj"
dim2 = ZobovTess(gal_rdz, gal, vol, zones, af, adj_filename)# initial class
dim2.make_voids(); dim2.make_zones()
#ej void 1
print dim2.zones_voids[1]
# above is the ID of zones that forms a watershed void
zones_void1 = dim2.zones_voids[1]
# print now the galaxies in each of these zones
print dim2.gals_zones[zones_void1[0]],dim2.gals_zones[zones_void1[1]],dim2.gals_zones[zones_void1[2]]
gals_in_void1 = np.concatenate(np.array(dim2.gals_zones)[zones_void1])
#plot void 1 in galaxies
V1 = dim2.gals_zobov[gals_in_void1]

#f = plt.figure(figsize=(8,6))
#ax = f.add_subplot(111,projection='3d')
#ax.plot(V1[:,0],V1[:,1],V1[:,2],'r+')
#plt.show()
