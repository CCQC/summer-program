import sys
sys.path.insert(0, '../../extra-files')

import numpy as np
from masses import get_charge, get_mass

class Molecule(object):

  def __init__(self, xyzstring, units="Angstrom"):
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
    self.assert_no_overlapping_atoms()

  def assert_no_overlapping_atoms(self):
    if not len(set(self.__str__().splitlines())) == len(list(self.__str__().splitlines())):
      raise Exception("Cannot construct Molecule with overlapping atoms.")

  def to_bohr(self):
    if self.units == "Angstrom":
      self.units  = "Bohr"
      self.geom  *= 1.889725989

  def to_angstrom(self):
    if self.units == "Bohr":
      self.units  = "Angstrom"
      self.geom  /= 1.889725989

  def copy(self):
    return Molecule(self.__str__(), self.units)

  def __len__(self):
    return self.natom

  def __iter__(self):
    for label, coords in zip(self.labels, self.geom):
      yield label, coords

  def __str__(self):
    fmt = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
    ret = "{:d}\n{:s}\n".format(self.natom, self.units)
    for label, coords in self.__iter__():
      ret += fmt.format(label, *coords)
    return ret

  def __add__(self, other):
    other = other.copy()
    if   self.units == "Angstrom" and other.units == "Bohr":
      other.to_angstrom()
    elif self.units == "Bohr"     and other.units == "Angstrom":
      other.to_bohr()
    other.natom   = self.natom   + other.natom
    other.labels  = self.labels  + other.labels
    other.masses  = self.masses  + other.masses
    other.charges = self.charges + other.charges
    other.geom    = np.concatenate((self.geom, other.geom))
    other.assert_no_overlapping_atoms()
    return other

    
    
      

if __name__ == "__main__":
  xyzstring = open("../../extra-files/molecule.xyz").read()
  mol = Molecule(xyzstring)
  print mol.units
  print mol.natom
  print mol.labels
  print mol.masses
  print mol.charges
  print mol.geom

  print(mol)

  mol.to_angstrom()
  print mol.units
  print mol.geom

  mol.to_bohr()
  print mol.units
  print mol.geom

  mol.to_angstrom()

  mol2 = Molecule(xyzstring)
  mol2.geom += 1

  print mol + mol2
  print mol + mol
