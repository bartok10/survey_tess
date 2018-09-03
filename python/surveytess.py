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
    self.zones_voids = []
    self.gals_zones = []
    self.Nzones = []
    self.volvoids = []
    self.gal_table = self.create_gal_table()
    self.adj_table = self.read_adj_ascii(adj_file)
    # self.void_table = self.create_void_table()

  def create_gal_table(self):
    """
      Method to save al ZOBOV galaxy information in tables
      """

    gal_table = Table()
    gal_table["X"] = self.gals_zobov[:,0]
    gal_table["Y"] = self.gals_zobov[:,1]
    gal_table["Z"] = self.gals_zobov[:,2]
    gal_table["vol_cell"] = self.vols_zobov
    gal_table['zone_id'] = self.zones_zobov.astype(int)
    gal_table['ID'] = range(0, len(gal_table))
    gal_table["RA"] = self.gals_orig[:,0]
    gal_table['DEC'] = self.gals_orig[:,1]
    gal_table['zgal'] = self.gals_orig[:,2]
    gal_table["RA_d"] = gal_table["RA"] * 180. / np.pi
    gal_table['DEC_d'] = gal_table["DEC"] * 180. / np.pi


    return gal_table

  def nearest_gal_neigh(self, p, coordinates='degree'):
      """
        :parameter:
        p : list
            [RA, DEC, z] with RA and DEC in degrees

        Find the closest galaxy to a point in comoving coordinates.
        Point "p" should be in (ra,dec) degrees coord and redshift as a float. Comoving distance
        is also allowed. By using closest neighbor algorithm with k-d tree
        optimization (in scipy module) a galaxy will be found.
        :return:
        The function return the ID of the neighbor galaxy.
       """
      if coordinates == 'degree':
          p_com = utils.deg2com(p)
      else:
          p_com = p
      pandgals = np.vstack([p_com, self.gals_zobov])
      # import pdb;pdb.set_trace()
      kdt = cKDTree(pandgals)
      p_neigh = kdt.query(pandgals,k=2)
      print("Distance to closest neighbour is {:.1f} Mpc".format(p_neigh[0][0,1]))
      ind = p_neigh[1][0,1] - 1  # -1 because p adds 1 to all indices
      if ind == -1:
        ind = p_neigh[1][0,0] - 1  # this is something weird of kdt.query() sometimes it flips the indices
      return ind

  def create_void_table(self):
    void_table = Table()
    void_table["void_gals"] = self.voids_zobov
    void_table["zones"] = self.zones_voids
    void_table["n_zones"] = self.Nzones
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

  def make_zones(self):
    """
    List of zones from ZOBOV with galaxy members.
    :return:
    table of zones with index of galaxies.
    """
    self.gals_zones = utils.cell2zones(self.zones_zobov)

  def make_voids(self):
    """
    List of void is built with the information given by ZOBOV.
    :return:
    table of voids with one or more zones formed by the Voronoi cells.
    """
    af=self.voids_zobov
    apf = []
    for i in range(len(af)):
        apf.append(af[i].split())
    vovp = []
    for i in range(1,len(apf)-1):
        vovp.append(utils.overlap(apf[i]))
    zone_IDs = []; nzones=[]
    for i in range(len(vovp)):
        if vovp[i] != -1:
            zone_IDs.append(vovp[i])
            nzones.append(len(vovp[i]))
        elif vovp[i] == -1:
            zone_IDs.append(i)
            nzones.append(0)
    cond = np.array([elem != 0 for elem in nzones])
    self.zones_voids = np.array(zone_IDs,dtype=object)[cond]
    self.Nzones = np.array(nzones)[cond]    


  def sbox_void(self,V):
      """
        Tessellation of small zone containing the void V and background galaxies to determine its properties.
        :param V: Watershed void with the position of the galaxy members.
        :return: Voronoi object 'vor' with all the information to be used in the voronoi_properties class.
        cv is a list with the IDs of the galaxies in the Voronoi object 'vor'.
      """
      xmax = max(V[:, 0]);ymax = max(V[:, 1]);zmax = max(V[:, 2])
      xmin = min(V[:, 0]);ymin = min(V[:, 1]);zmin = min(V[:, 2])
      bbx = (self.gals_zobov[:, 0] > xmin - abs(0.2 * xmin)) & (self.gals_zobov[:, 0] < xmax + abs(0.2 * xmax)) & (
                  self.gals_zobov[:, 2] > zmin - abs(0.2 * zmin)) & (self.gals_zobov[:, 2] < zmax + abs(0.2 * zmax)) & (
                        self.gals_zobov[:, 1] > ymin - abs(0.2 * ymin)) & (self.gals_zobov[:, 1] < ymax + abs(0.2 * ymax))
      sbox = np.array([self.gals_zobov[bbx, 0], self.gals_zobov[bbx, 1], self.gals_zobov[bbx, 2]]).T
      #vor = Voronoi(sbox)
      #cv = [] #cv contains the galaxy points IDs inside a void in "vor" space
      #for i in range(len(V)):
      #    l = len(np.where(vor.points[:, 0] == V[i][0])[0])
      #    if l != 0: cv.append(np.where(V[i][0] == vor.points[:, 0])[0][0])
      #########################
      #### select and label all the galaxies inside the void
      return sbox # voronoi scipy class

class voronoi_properties(): #voronoi properties of a void
    """
        Class to determine the Voronoi properties of a single galaxy:

        :parameters:
            vor : Voronoi object from the Voronoi class in scipy.spatial module

    """
    def __init__(self,voro):
        self.vor_tess = voro
        self.id = 0
        self.vol = 0.
        self.region = 0
        self.vertices = []
        self.faces = []
        self.adjs = []

    def vor_gal_id(self,gid):
        self.id = gid
    #def cell_vertices(self):
        self.region = self.vor_tess.point_region[self.id] #region ID for galaxy ID
        v_pol = self.vor_tess.regions[self.region] # vertices IDs for region ID
        vp = self.vor_tess.vertices[v_pol] #vertices positions
        self.vertices = vp

    #def cell_faces(self):
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


    #def cell_vol(self):
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




