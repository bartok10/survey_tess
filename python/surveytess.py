# Class for defining a Tessellated Survey of galaxies
import numpy as np
from scipy.spatial import Voronoi
from astropy.table import Table
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

  def create_void_table(self):
    void_table = Table()
    return void_table()

  def get_voronoi_cel_from_id(self, gal_id, plot=False, **kwargs):
    cond = self.adj_table['gal_id'] == gal_id
    assert np.sum(cond) == 1, "Something is wrong with the GAL ID in ADJ table"
    ids_adj = self.adj_table[cond]['ids_adj'].data[0]
    if ids_adj is None:
      return
    points_ids = [gal_id]
    for kk in ids_adj:
      points_ids += [kk]

    points = []
    for id in points_ids:
      x = self.gal_table[id]["X"]
      y = self.gal_table[id]["Y"]
      z = self.gal_table[id]["Z"]
      points += [[x,y,z]]
    # import pdb; pdb.set_trace()
    vor = Voronoi(points, **kwargs)  # Scipy Voro
    if plot:
      utils.plot_voronoi(vor)

    return vor



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


