# coding: utf-8
import numpy as np
from surveytess import  ZobovTess
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors

voids_id = np.array(dim2.gal_table['zone_id'])
void_id = dim2.gal_table['zone_id'][6]
gal2Vor = np.where(voids_id.astype(int) == void_id)[0]
f = plt.figure(figsize=(8,6))
ax = f.add_subplot(111,projection='3d')
for i in range(len(gal2Vor)):
    a = dim2.get_voronoi_cell_from_id(gal2Vor[i])
    gal_pos = a.points[0]
    cFaces = dim2.vorocell(a)
    ax.plot([gal_pos[0]], [gal_pos[1]], [gal_pos[2]], 'ro',markersize=4.)
    for i in range(len(cFaces)):
        poli = Poly3DCollection([a.vertices[cFaces[i]]])
        poli.set_edgecolor('k')
        poli.set_facecolor('grey')
        poli.set_alpha(alpha=0.1)
        ax.add_collection3d(poli)
