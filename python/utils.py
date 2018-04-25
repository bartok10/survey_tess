import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors


def gals2zones(zones, gal, vol):
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
    """
    :param p: (RA,DEC,z)
     RA and Dec in degrees
    :return:
    return (X,Y,Z) in comoving Mpc
    """
    ra_r = np.radians(p[0])
    dec_r = np.radians(p[1])
    z = p[2]
    dc_z = cosmo.comoving_distance(z).to('Mpc').value
    xp = dc_z * np.cos(ra_r) * np.cos(dec_r)
    yp = dc_z * np.sin(ra_r) * np.cos(dec_r)
    zp = dc_z * np.sin(dec_r)
    return np.array([xp,yp,zp])

