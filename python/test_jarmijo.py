import numpy as np
from surveytess import  ZobovTess,voronoi_properties
from scipy.spatial import Voronoi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
##################### solo en windows  ######################
#import sys
#sys.path.insert(0, 'C:/Users/Joaquin/survey_tess/python/')
#dir1 = 'C:/Users/Joaquin/Documents/test_st/'
#gal_rdz = np.loadtxt(dir1+'subsample.txt')
#gal = np.loadtxt(dir1+'gal_pos_3d_com.txt')
#vol = np.loadtxt(dir1+'cat_zobov.ascii.vol')
#zones = np.loadtxt(dir1+'cat_zobov.zone')
#f = open(dir1+'cat_zobov.void','r')
#af  = f.read().splitlines()
#f.close()
#adj_filename = dir1+"cat_zobov.ascii.adj"

##############################################################
#load files
# dir1 = '/home/jarmijo/CHUVIS/'
#dir1 =  '/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/'
dir1 = '/home/ntejos/python/survey_tess/CHUVIS/'

gal = np.loadtxt(dir1+'gal_XYZ.dat')
gal_rdz = np.loadtxt(dir1+'subsample.dat')
vol = np.loadtxt(dir1+'cat_zobov.ascii.vol',skiprows=1)
zones = np.loadtxt(dir1+'cat_zobov.zone',skiprows=1)
f = open(dir1+'cat_zobov.void','r')
af  = f.read().splitlines()
f.close()
adj_filename = dir1+"cat_zobov.ascii.adj"
################################################################

dim2 = ZobovTess(gal_rdz, gal, vol, zones, af, adj_filename)# initial class

dim2.make_voids() #method to build the watershed voids using the ZOBOV data

V = dim2.voids_zobov[236] #void de prueba

Vbox = dim2.sbox_void(V)

voro_void = Voronoi(Vbox)

cv = [] #cv contains the galaxy points IDs inside a void in "vor" space
for i in range(len(V)):
  l = len(np.where(voro_void.points[:, 0] == V[i][0])[0])
  if l != 0: cv.append(np.where(V[i][0] == voro_void.points[:, 0])[0][0])

celda = voronoi_properties(voro_void)
#celda.vor_gal_id(cv[2])
#
vol_void = []
for i in range(len(cv)):
    celda.vor_gal_id(cv[i])
    vol_void.append(celda.vol)
vol_V = sum(vol_void)
r_void = np.cbrt(3./4./np.pi*vol_V)
#
f = plt.figure(figsize=(8,6))
ax = f.add_subplot(111,projection='3d')
points = celda.vor_tess.points[cv]
bg = celda.vor_tess.points
ax.plot(points[:,0],points[:,1],points[:,2],'r+')
ax.plot(bg[:,0],bg[:,1],bg[:,2],'k.',ms=1.)
for i in range(len(cv)):
    celda.vor_gal_id(cv[i])
    faces = celda.faces
    for j in range(len(faces)):
        poli = Poly3DCollection([faces[j]])
        poli.set_edgecolor('black')
        poli.set_facecolor('red')
        ax.add_collection(poli)
plt.show()

def test_nearest_gal_neigh():
    for i in range(len(dim2.gals_zobov)):
        assert i == dim2.nearest_gal_neigh(dim2.gals_zobov[i], coordinates='com')
