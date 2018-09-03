import numpy as np
from surveytess import  ZobovTess

"""Module for testing the methods of ZobovTess"""


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

# def test_nearest_gal_neigh():
if 1:
    p = [149.21800,41.25610, 0.04377]
    dim2.nearest_gal_neigh(p, coordinates='degree')

    stop
    for i in range(len(dim2.gals_zobov)):
        assert i == dim2.nearest_gal_neigh(dim2.gals_zobov[i], coordinates='com')
