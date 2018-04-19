import numpy as np
from surveytess import  ZobovTess,voronoi_properties
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors

#dir1 = '/home/jarmijo/test_st/'
dir1 = '/home/joaquin/test_st/'
#dir1 =  '/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/'
gal = np.loadtxt(dir1+'gal_pos.txt')
vol = np.loadtxt(dir1+'cat_zobov.ascii.vol')
zones = np.loadtxt(dir1+'cat_zobov.zone')

f = open(dir1+'cat_zobov.void','r')
af  = f.read().splitlines()
f.close()

adj_filename = dir1+"cat_zobov.ascii.adj"

dim2 = ZobovTess([], gal, vol, zones, af, adj_filename)# initial class

V = dim2.make_voids() #method to build the watershed voids using the ZOBOV data
VV = V[7] #choose a void


