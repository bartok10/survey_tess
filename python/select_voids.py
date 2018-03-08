# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_void_from_id(id, stvoids, ax):
    
    xp = stvoids[id][:,0]
    yp = stvoids[id][:,1]
    zp = stvoids[id][:,2]
    # ax.plot(xe,ye,ze,'g+',markersize=0.2)
    # ax.plot(xc1,yc1,zc1,'b+',markersize=0.2)
    # ax.plot(xc2,yc2,zc2,'r+',markersize=0.2)
    ax.plot(xp,yp,zp,'o',markersize=5)
    ax.set_xlabel('X',fontsize=20)
    ax.set_ylabel('Y',fontsize=20)
    ax.set_zlabel('Z',fontsize=20)
    
    
def plot_all_gal(gal, ax):
    xp = gal[:,0]
    yp = gal[:,1]
    zp = gal[:,2]
    # ax.plot(xe,ye,ze,'g+',markersize=0.2)
    # ax.plot(xc1,yc1,zc1,'b+',markersize=0.2)
    # ax.plot(xc2,yc2,zc2,'r+',markersize=0.2)
    ax.plot(xp,yp,zp,'k+',markersize=2)
    ax.set_xlabel('X',fontsize=20)
    ax.set_ylabel('Y',fontsize=20)
    ax.set_zlabel('Z',fontsize=20)
    # plt.show()


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
apf = []
for i in range(len(af)):
    apf.append(af[i].split())
stvoids = []
stvol = []
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
vovp = []
for i in range(1,len(apf)):
    vovp.append(overlap(apf[i]))
    
# %load 100
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


# ploting
for ii in range(len(stvol_new[:10])):
    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(111,projection='3d')
    plot_all_gal(gal, ax)
    plot_void_from_id(ii, stvoids_new, ax)
    plt.show()
