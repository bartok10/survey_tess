import numpy as np
from surveytess import  ZobovTess
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors

dir1 = '/home/joaquin/test_st/'
#dir1 =  '/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/'
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

dim2 = ZobovTess([], gal, vol, zones, af, adj_filename)

a = dim2.get_voronoi_cel_from_id(0)#recibe ID galaxia
vert = np.array(a.vertices)
ridges = np.array(a.ridge_vertices) #caras de un poliedro

faces=[] #caras sin -1
for ll in ridges:
    ct=0
    for c in ll:
        if c==-1: break
        ct +=1
    if ct == len(ll):
        faces.append(ll)

#script para plottear celdas Voronoi (utils.py)
f = plt.figure(figsize=(8,6))
ax = f.add_subplot(111,projection='3d')
gal_pos = a.points[0]
adj_pos = a.points[1:]
xp = [pos[0] for pos in adj_pos]
yp = [pos[1] for pos in adj_pos]
zp = [pos[2] for pos in adj_pos]
scale=1.0
ax.plot(xp,yp,zp,'ob',markersize=2*scale)
ax.set_xlabel('X',fontsize=20)
ax.set_ylabel('Y',fontsize=20)
ax.set_zlabel('Z',fontsize=20)
ax.plot([gal_pos[0]], [gal_pos[1]], [gal_pos[2]], 'ro',markersize=4*scale)
vtx = [pos[0] for pos in a.vertices]
vty = [pos[1] for pos in a.vertices]
vtz = [pos[2] for pos in a.vertices]
ax.plot(vtx, vty, vtz, '^g',markersize=4*scale)
for i in range(len(faces)):
    poli = Poly3DCollection([a.vertices[faces[i]]])
    poli.set_edgecolor('k')
    poli.set_facecolor('gray')
    poli.set_alpha(alpha=0.1)
    ax.add_collection3d(poli)