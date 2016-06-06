#!/usr/bin/python 

import sys 
import numpy as np

sys.path.insert(0,'../../extra-files')
import masses as M

File = open('molecule.xyz', 'r')
File = File.read()
class Molecule(object):
   def __init__(self, xyz, units):
      self.units = units
      self.natom = int(xyz[0])
      self.xyz = []
      temp = xyz.splitlines()
      temp = temp[2:self.natom+2]
      for i in temp:
         self.xyz.append(i.split())
      self.labels = [i[0] for i in self.xyz]
      self.masses = [float(M.get_mass(self.labels[i])) for i in range(self.natom)]      
      self.charges = [M.get_charge(self.labels[i]) for i in range(self.natom)]
      geom = [[float(i[1]),float(i[2]),float(i[3])] for i in self.xyz]
      self.geom = np.matrix(geom) #This stores the list of lists as a matrix
   
   def to_bohr(self):
      if self.units == 'Angstrom':
         self.geom /= 0.529177208 
         self.units = 'Bohr'
      else:
         print('Units are already in Bohr.')

   def to_angstrom(self):
      if self.units == 'Bohr':
         self.geom *= 0.529177208
         self.units = 'Angstrom'

   def xyz_string(self):
      string = '%d \n%s \n' %(self.natom, self.units)
      for i in range(self.natom):
         string += '%s%20.10f%20.10f%20.10f\n' %(self.labels[i],self.geom[i,0],self.geom[i,1],self.geom[i,2])
      return string
   
   def copy(self):
      copy = Molecule(self.xyz_string(),self.units)
      return copy   

mol = Molecule(File, 'Angstrom')

print(mol.xyz)
print(mol.labels)
print(mol.masses)
print(mol.charges)
print(mol.geom)
print(mol.xyz_string())

mol2 = mol.copy()
print mol2.masses

