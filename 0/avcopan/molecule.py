import sys
sys.path.insert(0, '../../extra-files')

import numpy as np
from masses import get_charge, get_mass

class Molecule(object):

  def __init__(self, xyzstring, units="angstrom"):
    labels = []
    geom   = []
    try:
      lines  = xyzstring.strip().splitlines()
      natom  = int(lines[0])
      for line in lines[2:]:
        l, x, y, z = line.split()
        labels.append(l.upper())
        geom.append([float(x), float(y), float(z)])
    except:
      raise Exception("Invalid .xyz string\n'''\n{:s}\n'''\npassed to Molecule constructor.".format(xyzstring))
    self.units   = units
    self.natom   = natom    
    self.labels  = labels
    self.masses  = [get_mass(label)   for label in labels]
    self.charges = [get_charge(label) for label in labels]
    self.geom    = np.array(geom)

  def to_bohr(self):
    if self.units == "angstrom":
      self.units  = "bohr"
      self.geom  *= 1.889725989

  def to_angstrom(self):
    if self.units == "bohr":
      self.units  = "angstrom"
      self.geom  /= 1.889725989

  def copy(self):
    return Molecule(self.xyz_string(), self.units)

  def xyz_string(self):
    ret = "{:d}\nunits {:s}\n".format(self.natom, self.units)
    fmt = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
    for label, coords in zip(self.labels, self.geom):
      ret += fmt.format(label, *coords)
    return ret


if __name__ == "__main__":
  xyzstring = open("../../extra-files/molecule.xyz").read()
  mol = Molecule(xyzstring)
  print mol.units
  print mol.natom
  print mol.labels
  print mol.masses
  print mol.charges
  print mol.geom

  mol.to_angstrom()
  print mol.xyz_string()

  mol.to_bohr()
  print mol.xyz_string()

  mol2 = mol.copy()
  print mol.xyz_string()
