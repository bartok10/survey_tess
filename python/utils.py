import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
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
    voids = []
    if len(void) > 3.:
        j = 1
        cont = 3
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
    #plt.show()

def deg2com(p): #degree + redshift 3D coordinates
    ra = p[0]; dec = p[1]; z = p[2]
    dc_z = cosmo.comoving_distance(z)
    xp = dc_z * np.cos(ra) * np.sin(dec)
    yp = dc_z * np.sin(ra) * np.sin(dec)
    zp = dc_z * np.cos(dec)
    return np.array([xp,yp,zp])
