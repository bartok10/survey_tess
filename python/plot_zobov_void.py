# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
from scipy.spatial import Voronoi
#
gal = np.loadtxt('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/gal_pos.txt')
vol = np.loadtxt('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/cat_zobov.ascii.vol')
zones = np.loadtxt('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/cat_zobov.zone')
v_mem = []
v_vols = []
#
for i in range(int(max(zones))+1):
    v_mem.append(gal[zones==i])
    v_vols.append(vol[zones==i])
def overlap(void):
    #v_0 = int(void[0])
    if len(void) > 3.:
        j = 1
        cont = 3
        voids = []
        while int(void[j]) != 0:
            v_j = np.zeros(int(void[j]))
            for i in range(len(v_j)):
                voids.append(int(void[cont+i]))
            cont = cont + len(v_j) + 2
            j = j + len(v_j) + 2
    elif len(void) == 3.:
        voids = -1
    return voids
f = open('/media/ntejos/disk1/projects/COS-Web/voronoi/test_st/cat_zobov.void','r')
af  = f.read().splitlines()
f.close()
apf = []
for i in range(len(af)):
    apf.append(af[i].split())
stvoids = []
stvol = []
vovp = []
for i in range(1,len(apf)):
    vovp.append(overlap(apf[i]))
for i in range(len(vovp)):
    if vovp[i] != -1 :
        tmpv=[];tmpvl=[]
        for j in range(len(vovp[i])):
            tmpv.append(v_mem[vovp[i][j]])
            tmpvl.append(v_vols[vovp[i][j]])
        tmpa = np.concatenate(tmpv,axis=0)
        tmpb = np.concatenate(tmpvl,axis=0)
        stvoids.append(np.concatenate([v_mem[i],tmpa]))
        stvol.append(np.concatenate([v_vols[i],tmpb]))
    elif vovp[i] == -1 :
        stvoids.append(v_mem[i])
        stvol.append(v_vols[i])
stvoids = np.array(stvoids)
stvol = np.array(stvol)
cond = np.array([elem != -1 for elem in vovp])

stvoids_new = stvoids[cond]
stvol_new = stvol[cond]

V =stvoids_new[19]
xmax = max(V[:,0]); ymax = max(V[:,1]); zmax = max(V[:,2])
xmin = min(V[:,0]); ymin = min(V[:,1]); zmin = min(V[:,2])
bbx = (gal[:,0]>xmin-abs(0.1*xmin))&(gal[:,0]<xmax+abs(0.1*xmax))&(gal[:,2]>zmin-abs(0.1*zmin))&(gal[:,2]<zmax+abs(0.1*zmax))&(gal[:,1]>ymin-abs(0.1*ymin))&(gal[:,1]<ymax+abs(0.1*ymax))
sbox = np.array([gal[bbx,0],gal[bbx,1],gal[bbx,2]]).T
vor = Voronoi(sbox)
cv = []
for i in range(len(V)):
    l = len(np.where(vor.points[:,0]==V[i][0])[0])
    if l != 0: cv.append(np.where(V[i][0]==vor.points[:,0])[0][0])
   

#
f = plt.figure(figsize=(8,6))
ax = f.add_subplot(111,projection='3d')
ax.plot(V[:,0],V[:,1],V[:,2],'g*')
ax.plot(gal[bbx,0],gal[bbx,1],gal[bbx,2],'k.',markersize=0.5)
plt.show()
####################################

####################################

#
p = cv[0]
pc = vor.points[p]
pol_id = vor.point_region[p]
v_pol = vor.regions[pol_id]
vp = vor.vertices[v_pol]
f = plt.figure(figsize=(8,6))
ax = f.add_subplot(111,projection='3d')
ax.plot([pc[0]],pc[1],pc[2],'r*')
ax.plot(vp[:,0],vp[:,1],vp[:,2],'g+')
plt.show()
#
ridges_id = []
ridges = vor.ridge_vertices
for i in range(len(ridges)):
    for j in range(len(cv)):
        if vor.ridge_points[i,0] == cv[j]:
            ridges_id.append(i); break
#         
f = plt.figure(figsize=(8,6))
ax = f.add_subplot(111,projection='3d')
pc = vor.points[cv[0]]
for i in range(len(ridges_id)):
    face_i = ridges_id[i]
    vface_i = vor.ridge_vertices[face_i]
    poli = Poly3DCollection([vor.vertices[vface_i]])
    poli.set_edgecolor('white')
    poli.set_facecolor('grey')
    poli.set_alpha(alpha=0.1)
    ax.add_collection3d(poli)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.plot(gal[bbx,0],gal[bbx,1],gal[bbx,2],'k.',markersize=0.5)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_zlim(zmin,zmax)
plt.show()
              
       
#get_ipython().magic(u'save plot_void_dbg.py ')
