# A Molecule class that takes a zmatrix as its argument.  Internally, it builds a
# Cartesian representation of the geometry (mass-centered, with x, y, z aligned to
# the inertial axes as I_x < I_y < I_z) and uses this to compute angles, direction
# vectors, etc.  One of the methods, self.s_vector, returns s_vectors for distances,
# angles, and torsions.  For example, self.s_vector(a, (a, b)) returns s_a(r_ab),
# self.s_vector(a, (a, b, c)) returns s_a(phi_abc), and self.s_vector(a, (a, b, c, d))
# returns s_a(tau_abcd).  It also has the method self.compute_g_matrix() which uses
# the s vectors to build the G matrix and return it.  Depends on the module "masses"
# which simply contains all of the isotopic masses and returns them by atomic symbol
# through the function "get_mass()".
import numpy as np
import scipy.linalg as la
from masses import get_mass

class Molecule(object):

  def __init__(self, zmatrix):
    labels = []
    zmat_keys = []
    zmat = {}
    for atom1, line in enumerate(zmatrix.strip().splitlines()):
      items = line.split()
      labels.append( items.pop(0) )
      key = (atom1,)
      for atom2, value in zip(items[0::2], items[1::2]):
        key += (int(atom2) - 1,) # -1 because we're 0-indexing
        zmat_keys.append(key)
        value = float(value) if len(key) < 3 else float(value) * 2 * np.pi / 360. # convert to radians if it's an angle
        zmat.update({key: float(value)})
    self.natoms    = len(labels)
    self.nintcos   = len(zmat_keys)
    self.labels    = labels
    self.masses    = [get_mass(label) for label in labels]
    self.zmat_keys = zmat_keys
    self.zmat      = zmat
    self.geom      = np.zeros((self.natoms, 3))
    self.fill_geom_from_zmat()

  def fill_geom_from_zmat(self):
    keys = self.zmat_keys[:]
    ex, ey, ez = np.identity(3)
    fabc = np.pi       # overwritten after 0, 1
    tabcd = np.pi / 2  # overwritten after 0, 1, 2
    ecb = ex           # overwritten after 0, 1
    ecd = ez           # overwritten after 0, 1, 2
    try:
      for atom in range(1, self.natoms):
        a, b = keys.pop(0)
        rab = self.zmat[(a, b)]
        if atom >= 2:
          a, b, c = keys.pop(0)
          fabc = self.zmat[(a, b, c)]
          ecb = self.unit_direction_vector(c, b)
        if atom >= 3:
          a, b, c, d = keys.pop(0)
          ecd = self.unit_direction_vector(c, d)
          tabcd = self.zmat[(a, b, c, d)]
        ndcb = np.cross(ecb, ecd) / la.norm( np.cross(ecb, ecd) )
        nper = np.cross(ndcb, ecb)
        self.geom[a] = self.geom[b] \
                     - rab * np.cos(fabc) * ecb \
                     + rab * np.sin(fabc) * np.sin(tabcd) * ndcb \
                     + rab * np.sin(fabc) * np.cos(tabcd) * nper 
      self.geom = self.geom - self.compute_center_of_mass()
      self.moments, R = self.diagonalize_inertia_tensor()
      self.geom = np.dot(self.geom, R)
    except:
      raise Exception("Invalid zmatrix passed to Molecule constructor.")

  def compute_center_of_mass(self):
    M = sum(self.masses)
    return sum(m*r for m, r in zip(self.masses, self.geom)) / M

  def diagonalize_inertia_tensor(self):
    I = sum( m*(np.dot(r, r)*np.eye(3) - np.outer(r, r)) for m, r in zip(self.masses, self.geom) )
    return la.eigh(I)

  def distance(self, a, b):
    return la.norm( self.direction_vector(a, b) )

  def angle(self, a, b, c):
    e = self.unit_direction_vector
    return np.arccos(np.dot(e(b, a), e(b, c)))

  def torsion(self, a, b, c, d):
    e = self.unit_direction_vector
    f = self.angle
    return np.arccos(np.dot(np.cross(e(a, b), e(b, c)), np.cross(e(b, c), e(c, d))) / (np.sin(f(a, b, c)) * np.sin(f(b, c, d))))

  def direction_vector(self, a, b):
    return self.geom[b] - self.geom[a]

  def unit_direction_vector(self, a, b):
    return self.direction_vector(a, b) / self.distance(a, b)

  def unit_normal_vector(self, a, b, c):
    e = self.unit_direction_vector
    f = self.angle
    return np.cross(e(b, a), e(b, c)) / np.sin(f(a, b, c))

  def s_vector(self, atom, zmat_key):
    r = self.distance
    f = self.angle
    e = self.unit_direction_vector
    n = self.unit_normal_vector
    if not atom in zmat_key:
      return np.zeros(3)
    if len(zmat_key) == 2:
      a, b = zmat_key
      return e(b, a) if atom is a else e(a, b)
    if len(zmat_key) == 3:
      a, b, c = zmat_key
      sa = ( e(b, a)*np.cos(f(a, b, c)) - e(b, c) ) / ( r(b, a)*np.sin(f(a, b, c)) )
      sc = ( e(b, c)*np.cos(f(a, b, c)) - e(b, a) ) / ( r(b, c)*np.sin(f(a, b, c)) )
      return sa if atom is a else ( sc if atom is c else -sa -sc)
    if len(zmat_key) == 4:
      a, b, c, d = zmat_key
      if atom in (c, d): a, b, c, d = d, c, b, a
      if atom is a:
        return n(a, b, c) / ( r(a, b)*np.sin(f(a, b, c)) )
      elif atom is b:
        return n(c, b, a)*(r(b, c) - r(a, b)*np.cos(f(a, b, c)))/(r(a, b)*r(b, c)*np.sin(f(a, b, c))) \
             + n(b, c, d)* np.cos(f(b, c, d))/(r(b, c)*np.sin(f(b, c, d)))

  def compute_g_matrix(self):
    s = self.s_vector
    m = self.masses
    g = np.zeros((self.nintcos, self.nintcos))
    for k1, key1 in enumerate(self.zmat_keys):
      for k2, key2 in enumerate(self.zmat_keys):
        g[k1, k2] = sum( np.dot(s(a, key1), s(a, key2)) / m[a] for a in set(key1) & set(key2) )
    return g
        
