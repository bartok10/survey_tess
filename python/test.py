# coding: utf-8
import numpy as np
from surveytess import ZobovTess

gal = np.loadtxt('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/gal_pos.txt')
vol = np.loadtxt('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/cat_zobov.ascii.vol')
zones = np.loadtxt('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/cat_zobov.zone')
v_mem = []
v_vols = []
for i in range(int(max(zones))+1):
    v_mem.append(gal[zones==i])
    v_vols.append(vol[zones==i])
f = open('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/cat_zobov.void','r')
af  = f.read().splitlines()
f.close()

adj_filename = "/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/cat_zobov.ascii.adj"

dim2 = ZobovTess([], gal, vol, zones, af, adj_filename)
