import numpy as np
from surveytess import  ZobovTess,voronoi_properties
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
##################### solo en windows  ######################
import sys
sys.path.insert(0, 'C:/Users/Joaquin/survey_tess/python/')
dir1 = 'C:/Users/Joaquin/Documents/test_st/'
gal_rdz = np.loadtxt(dir1+'subsample.txt')
gal = np.loadtxt(dir1+'gal_pos_3d_com.txt')
vol = np.loadtxt(dir1+'cat_zobov.ascii.vol')
zones = np.loadtxt(dir1+'cat_zobov.zone')
f = open(dir1+'cat_zobov.void','r')
af  = f.read().splitlines()
f.close()
adj_filename = dir1+"cat_zobov.ascii.adj"

##############################################################
#load files
#dir1 = '/home/joaquin/test_st/'
#dir1 =  '/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/'
#gal = np.loadtxt(dir1+'gal_pos.txt')
#vol = np.loadtxt(dir1+'cat_zobov.ascii.vol')
#zones = np.loadtxt(dir1+'cat_zobov.zone')
#f = open(dir1+'cat_zobov.void','r')
#af  = f.read().splitlines()
#f.close()
#adj_filename = dir1+"cat_zobov.ascii.adj"
################################################################

dim2 = ZobovTess(gal_rdz, gal, vol, zones, af, adj_filename)# initial class

V = dim2.make_voids() #method to build the watershed voids using the ZOBOV data

