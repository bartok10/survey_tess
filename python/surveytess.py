# Class for defining a Tessellated Survey of galaxies
import numpy as np
from scipy.spatial import Voronoi
from astropy.table import Table
from scipy.spatial import cKDTree
import utils

class ZobovTess():
  """Class for handling a Voronoi tesselation of volume from a galaxy survey
  using the ZOBOV alrogithm

  Parameters
  ----------
  gals_orig :

  gals_zobov : numpy array
    List of all galaxy positions in co-moving coordinates
  vols_zobov : numpy array
    List of all Voronoi volumes of each galaxy in `gal` as given by ZOBOV
  zones_zobov : numpy array
    List of IDs of "Void Zones" as defined by ZOBOV for each galaxy in `gal`
  voids_zobov : numpy array
    List of individual "Voids" as defined by ZOBOV for each `zone`

  """

  def __init__(self, gals_orig, gals_zobov, vols_zobov, zones_zobov, voids_zobov, adj_file):
    self.gals_orig = gals_orig
    self.gals_zobov = gals_zobov
    self.vols_zobov = vols_zobov
    self.zones_zobov = zones_zobov
    self.voids_zobov = voids_zobov
    self.gal_table = self.create_gal_table()
    self.adj_table = self.read_adj_ascii(adj_file)

  def create_gal_table(self):
    gal_table = Table()
    gal_table["X"] = self.gals_zobov[:,0]
    gal_table["Y"] = self.gals_zobov[:,1]
    gal_table["Z"] = self.gals_zobov[:,2]
    gal_table["vol_cell"] = self.vols_zobov
    gal_table['zone_id'] = self.zones_zobov
    gal_table['ID'] = range(0,len(gal_table))
    # gal_table["RA"] = self.gals_orig["RA"]
    # gal_table['DEC'] = self.gals_orig["DEC"]
    # gal_table['zgal'] = self.gals_orig["z"]
    return gal_table

  def nearest_gal_neigh(self,p,coordinates='degree'):
      if coordinates == 'degree': p_com = utils.deg2com(p)
      else: p_com = p
      pandgals = np.vstack([p_com,self.gals_zobov])
      kdt = cKDTree(pandgals)
      p_neigh = kdt.query(pandgals,k=2)[1]
      return (p_neigh[0,1] - 1)

  def create_void_table(self):
    void_table = Table()
    return void_table

  def read_adj_ascii(self, ascii_adj):
    f = open(ascii_adj, 'r')
    lines = f.readlines()
    n_gals = int(lines[0][:-1])
    id_gal = []
    n_adj = []
    ids_adj = []
    for line in lines[n_gals+1:]:
      id_gal_aux = int(line.split(":")[0].split(" ")[0])
      n_adj_aux = int(line.split(":")[0].split(" ")[1])
      if n_adj_aux > 0:
        ids_adj_aux = line.split(":")[1].split(" ")[1:-1]
        ids_adj_aux = [int(elem) for elem in ids_adj_aux]
      else:  # this could happen probably ZOBOV's fault
        ids_adj_aux = None
      id_gal += [id_gal_aux]
      n_adj += [n_adj_aux]
      ids_adj += [ids_adj_aux]
    tab = Table()
    tab['gal_id'] = id_gal
    tab['n_adj'] = n_adj
    tab['ids_adj'] = ids_adj
    f.close()
    return tab

  def make_voids(self):
    cmem,cvols = utils.cell2zones(self.zones_zobov,self.gals_zobov,self.vols_zobov)
    af=self.voids_zobov
    apf = []
    for i in range(len(af)):
        apf.append(af[i].split())
    stvoids = []
    stvol = []
    vovp = []
    for i in range(1,len(apf)):
        vovp.append(utils.overlap(apf[i]))
    for i in range(len(vovp)):
        if vovp[i] != -1:
            tmpv = []; tmpvl = []
            for j in range(len(vovp[i])):
                tmpv.append(cmem[vovp[i][j]])
                tmpvl.append(cvols[vovp[i][j]])
            tmpa = np.concatenate(tmpv, axis=0)
            tmpb = np.concatenate(tmpvl, axis=0)
            stvoids.append(np.concatenate([cmem[i], tmpa]))
            stvol.append(np.concatenate([cvols[i], tmpb]))
        elif vovp[i] == -1:
            stvoids.append(cmem[i])
            stvol.append(cvols[i])

    stvoids = np.array(stvoids)
    stvol = np.array(stvol)
    cond = np.array([elem != -1 for elem in vovp])

    stvoids_new = stvoids[cond]
    stvol_new = stvol[cond]
    return stvoids_new


  def vorocells(self,V): #recieve a void V
      ##### create subbox #####
      xmax = max(V[:, 0]);ymax = max(V[:, 1]);zmax = max(V[:, 2])
      xmin = min(V[:, 0]);ymin = min(V[:, 1]);zmin = min(V[:, 2])
      bbx = (self.gals_zobov[:, 0] > xmin - abs(0.1 * xmin)) & (self.gals_zobov[:, 0] < xmax + abs(0.1 * xmax)) & (
                  self.gals_zobov[:, 2] > zmin - abs(0.1 * zmin)) & (self.gals_zobov[:, 2] < zmax + abs(0.1 * zmax)) & (
                        self.gals_zobov[:, 1] > ymin - abs(0.1 * ymin)) & (self.gals_zobov[:, 1] < ymax + abs(0.1 * ymax))
      sbox = np.array([self.gals_zobov[bbx, 0], self.gals_zobov[bbx, 1], self.gals_zobov[bbx, 2]]).T
      vor = Voronoi(sbox)
      cv = [] #cv contains the galaxy points IDs inside a void in "vor" space
      for i in range(len(V)):
          l = len(np.where(vor.points[:, 0] == V[i][0])[0])
          if l != 0: cv.append(np.where(V[i][0] == vor.points[:, 0])[0][0])
      #########################
      #### select and label all the galaxies inside the void
      return vor,cv # voronoi scipy class

class voronoi_properties(): #voronoi properties of a void
    def __init__(self,vor):
        self.vor_tess = vor
        self.id = 0
        self.vol = 0.
        self.region = 0
        self.vertices = []
        self.faces = []
        self.adjs = []
    def galaxy_id(self,gid):
        self.id = gid
    def cell_vertices(self):
        self.region = self.vor_tess.point_region[self.id] #region ID for galaxy ID
        v_pol = self.vor_tess.regions[self.region] # vertices IDs for region ID
        vp = self.vor_tess.vertices[v_pol] #vertices positions
        self.vertices = vp

    def cell_faces(self):
        ridges_ids = []
        adjs = []
        ridges = self.vor_tess.ridge_vertices
        for i in range(len(ridges)):
            if self.vor_tess.ridge_points[i, 0] == self.id:
                ridges_ids.append(i)
                adjs.append(self.vor_tess.ridge_points[i,1])
        self.adjs = adjs
        vert_faces = []
        for i in range(len(ridges_ids)):
            face_i = ridges_ids[i]
            vface_i = self.vor_tess.ridge_vertices[face_i]
            vert_faces.append(self.vor_tess.vertices[vface_i])

        self.faces = vert_faces


    def cell_vol(self):
        p_adj = self.vor_tess.points[self.adjs]
        p = self.vor_tess.points[self.id]
        for f in range(len(self.faces)):
            h = 0.5 * np.linalg.norm((p - p_adj[f]))
            area = 0.  # triangle area
            for i in range(len(self.faces[f]) - 2):
                v0 = self.faces[f][0]
                v1 = self.faces[f][i + 1]
                v2 = self.faces[f][i + 2]
                d1 = np.linalg.norm((v0 - v1))
                d2 = np.linalg.norm((v0 - v2))
                d3 = np.linalg.norm((v2 - v1))
                s = 0.5 * (d1 + d2 + d3)  # semiperimeter
                area += np.sqrt(abs(s * (s - d1) * (s - d2) * (s - d3)))  # Heron's formula
            v = (1. / 3.) * area * h #tetrahedron volume (Heron)
            self.vol += v




