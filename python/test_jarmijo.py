import numpy as np
from surveytess import  ZobovTess
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors

# dir1 = '/home/joaquin/test_st/'
dir1 =  '/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/'
gal = np.loadtxt(dir1+'gal_pos.txt')
vol = np.loadtxt(dir1+'cat_zobov.ascii.vol')
zones = np.loadtxt(dir1+'cat_zobov.zone')
v_mem = []
v_vols = []
for i in range(int(max(zones))+1):
    v_mem.append(gal[zones==i])
    v_vols.append(vol[zones==i])
f = open(dir1+'cat_zobov.void','r')
af  = f.read().splitlines()
f.close()

adj_filename = dir1+"cat_zobov.ascii.adj"

dim2 = ZobovTess([], gal, vol, zones, af, adj_filename)# initial class

V = dim2.make_voids() #method to build the watershed voids using the ZOBOV data
VV = V[7] #choose a void
VVV = dim2.vorocell(VV) # vorocell recieves a void and return the vertices of each voroni cell in the void

f = plt.figure(figsize=(8,6))
ax = f.add_subplot(111,projection='3d')
for i in range(len(VVV)):
     ax.plot(VVV[i][:,0],VVV[i][:,1],VVV[i][:,2],'+')
