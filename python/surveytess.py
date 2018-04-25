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
    List of IDs for individual "Voids" as defined by ZOBOV for each `zone`

  """

  def __init__(self, gals_orig, gals_zobov, vols_zobov, zones_zobov, voids_zobov, adj_file):
    self.gals_orig = gals_orig
    self.gals_zobov = gals_zobov
    self.vols_zobov = vols_zobov
    self.zones_zobov = zones_zobov
    self.voids_zobov = voids_zobov
    # self.zones_voids = []  # gal_pos_in_zones ?

    self.gal_table = self.create_gal_table()
    self.adj_table = self.read_adj_ascii(adj_file)
    #self.void_table = self.create_void_table()

  def create_gal_table(self):
    """
      Method to save al ZOBOV galaxy information in tables
      """
    gal_table = Table()
    gal_table["X"] = self.gals_zobov[:,0]
    gal_table["Y"] = self.gals_zobov[:,1]
    gal_table["Z"] = self.gals_zobov[:,2]
    gal_table["vol_cell"] = self.vols_zobov
    gal_table['zone_id'] = self.zones_zobov
    gal_table['ID'] = range(0,len(gal_table))
    gal_table["RA"] = self.gals_orig[:,0]
    gal_table['DEC'] = self.gals_orig[:,1]
    gal_table['zgal'] = self.gals_orig[:,2]
    return gal_table

  def create_void_table(self):
      void_table = Table()
      # get relevant information
      gal_pos, zones, n_zones = self.make_voids()
      void_table['ID'] = range(0, len(zones))  # Unique void IDs
      void_table["gal_pos"] = gal_pos
      void_table["zones"] = zones
      void_table["n_zones"] = n_zones
      return void_table

  def create_zone_table(self):
      zone_table = Table()
      id_zones = range(0, len(np.unique(self.zones_zobov)))
      zone_table["ID"] = id_zones


  def nearest_gal_neigh(self,p,coordinates='degree'):
    """
    Find the closest galaxy to a point in comoving coordinates.
    Point "p" should be in (ra,dec, z) coordinates. Comoving position
    (X, Y, Z) is also allowed. By using closest neighbor algorithm
    with k-d tree optimization (in scipy module) a galaxy will be found.

    Parameters
    ----------
    p : np.array (N=3)
        Point representing (RA, Dec, z)
        or (X, Y, Z) if given in comoving coordinates
    coordinates : str
        Type of coordinate to perform the search
            if 'degree' -> assumes (RA, Dec, z)
            if 'comoving' -> assumes (X, Y, Z) in Mpc

    Returns
    -------
    id, dist : tuple
        The ID and distances of the closest galaxy to a given point p.
    """
    if coordinates == 'degree':
        p_com = utils.deg2com(p)
    elif coordinates == 'comoving':
        p_com = p
    else:
        raise ValueError('Not prepared for this type of `coordinates`.')
    p_and_gals = np.vstack([p_com, self.gals_zobov])
    kdt = cKDTree(p_and_gals)
    dist_neigh, p_neigh = kdt.query(p_and_gals, k=2)
    id_gal = (p_neigh[0,1] - 1)
    dist_gal = dist[0,1]
    return id_gal, dist_gal



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
    """
    List of void is built with the information given by ZOBOV.
    :return:
    Watershed voids with one or more zones formed by the Voronoi cells.
    """
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
    zone_IDs = []; nzones=[]
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
            zone_IDs.append(vovp[i])
            nzones.append(len(vovp[i]))
        elif vovp[i] == -1:
            stvoids.append(cmem[i])
            stvol.append(cvols[i])

    stvoids = np.array(stvoids)
    stvol = np.array(stvol)
    cond = np.array([elem != -1 for elem in vovp])

    stvoids_new = stvoids[cond]
    stvol_new = stvol[cond]
    # self.voids_zobov = stvoids_new
    # self.zones_voids = zone_IDs
    # self.Nzones = nzones
    return stvoids_new, zone_IDs, nzones



  def sbox_void(self, V):
      """
      Obtain a subsample of galaxy positions around ZOBOV void V

      Parameters
      ----------
      V :


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

class VoronoiCell():
    """
        Class that represents a Voronoi cell

        :parameters:
            voro : Voronoi object from the Voronoi class in scipy.spatial module
                 Has many galaxies, and we only want to characterize the Voronoi Cell associated
                 to 1 special galaxy
            gal_id : int
                The ID for the galaxy of interest within the voro object

        Atributes
        ---------

        Methods
        -------

    """
    def __init__(self, voro, gal_id, set_volume=True):
        self.voro = voro  # Voronoi object from scipy, has many galaxies beyond what is relevant for this class
        self.gal_id = gal_id  # ID of galaxy defining the individual Voronoi object of interest
        self.vor_id = self.voro.point_region[self.gal_id] # ID of the voronoi cell defined by galaxy gal_id
        self.gal_pos = self.voro.points[self.gal_id]  # position of the galaxy defining the voronoi cell
        self.vertices_ids = self.voro.regions[self.vor_id] # vertices IDs for voronoi cell vor_id ID
        self.vertices = self.voro.vertices[self.vertices_ids]  # Vertices positions of the voronoi cell of interest
        self.faces, self.adjs = self.cell_faces_and_adjacents()  # List of positions of vertices per face
        if set_volume:
            self.volume = self.get_volume()  # get the volume in Mpc^3
        else:
            self.volume = None

    def cell_faces_and_adjacents(self):
        """Get the vertices IDs per face"""
        ridges_ids = []
        adj_ids = []
        # first, get the IDs for ridges and adjs points
        ridges = self.voro.ridge_vertices  # list of "ridges" for voro object
        for i in range(len(ridges)):
            if self.voro.ridge_points[i, 0] == self.gal_id:
                ridges_ids.append(i)
                adj_ids.append(self.voro.ridge_points[i,1])
        # then, for each face of interest, get a list of vertices ids per face
        vert_faces = []
        for ids in ridges_ids:
            ids_vertices_in_face = self.voro.ridge_vertices[ids]  # indices of vertices of face in `ids`
            vert_faces.append(self.voro.vertices[ids_vertices_in_face])
        # return
        return vert_faces, adj_ids

    def get_volume(self):
        if self.volume is not None:
            return self.volume

        # if not defined already, calculate it with Heron's formula
        adj_points = self.voro.points[self.adjs]
        volume = 0.
        for f_ii, face in enumerate(self.faces):
            # h is the height between face f_ii and gal_pos
            h = 0.5 * np.linalg.norm((self.gal_pos - adj_points[f_ii]))
            area = 0.  # triangle area
            for f_jj in range(len(self.faces[f_ii]) - 2):
                # vertices of individual triangles
                v0 = self.faces[f_ii][0]
                v1 = self.faces[f_ii][f_jj + 1]
                v2 = self.faces[f_ii][f_jj + 2]
                # sides of those triangles
                d0 = np.linalg.norm((v0 - v1))
                d1 = np.linalg.norm((v0 - v2))
                d2 = np.linalg.norm((v2 - v1))
                # semi-perimeter and Heron formula
                s = 0.5 * (d0 + d1 + d2)  # semi-perimeter
                area += np.sqrt(abs(s * (s - d0) * (s - d1) * (s - d2)))  # Heron's formula
            v = (1. / 3.) * area * h  # tetrahedron volume (Heron)
            volume += v
        return volume

