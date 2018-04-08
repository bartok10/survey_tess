import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors


def cell2zones(zones,gal,vol):
    v_mem,v_vols = [],[]
    for i in range(int(max(zones))+1):
        v_mem.append(gal[zones==i])
        v_vols.append(vol[zones==i])
    return v_mem, v_vols

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


def plot_voronoi(vor, scale=1.):
    f = plt.figure(figsize=(8,6))
    ax = f.add_subplot(111,projection='3d')
    gal_pos = vor.points[0]
    adj_pos = vor.points[1:]
    xp = [pos[0] for pos in adj_pos]
    yp = [pos[1] for pos in adj_pos]
    zp = [pos[2] for pos in adj_pos]
    ax.plot(xp,yp,zp,'ob',markersize=2*scale)
    ax.set_xlabel('X',fontsize=20)
    ax.set_ylabel('Y',fontsize=20)
    ax.set_zlabel('Z',fontsize=20)
    ax.plot([gal_pos[0]], [gal_pos[1]], [gal_pos[2]], 'ro',markersize=4*scale)

    vtx = [pos[0] for pos in vor.vertices]
    vty = [pos[1] for pos in vor.vertices]
    vtz = [pos[2] for pos in vor.vertices]
    ax.plot(vtx, vty, vtz, '^g',markersize=4*scale)

    # identify the closest vertices to gal_pos
    #dist = np.sqrt((vtx-gal_pos[0])**2 + (vty-gal_pos[1])**2 + (vtz-gal_pos[2])**2)
    # euler's formula V - E + F = 2 (in 3D)
    #F = len(adj_pos)

    #poli = Poly3DCollection([vor.vertices])
    #poli.set_edgecolor('k')
    #ax.add_collection3d(poli)

    #import pdb; pdb.set_trace()
    #plt.show()