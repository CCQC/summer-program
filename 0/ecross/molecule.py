#!/usr/bin/python 

import sys 
import numpy as np


sys.path.insert(0,'../../extra-files')
import masses as M

class Molecule(object):
   def __init__(self,molecule, units):
      File = open(molecule, 'r')
      a = File.read()
      self.units = units
      self.natom = int(a[0])
      self.xyz = []
      temp = a.splitlines()
      temp = temp[2:self.natom+2]
      for i in temp:
         self.xyz.append(i.split())
      self.labels = [i[0] for i in self.xyz]
      self.masses = [float(M.get_mass(self.labels[i])) for i in range(self.natom)]      
      self.charges = [M.get_charge(self.labels[i]) for i in range(self.natom)]
      geom = [[float(i[1]),float(i[2]),float(i[3])] for i in self.xyz]
      self.geom = np.matrix(geom) #This stores the list of lists as a matrix
      File.close() 



   def to_bohr(self):   
#      if self.units != 'Bohr' or 'Angstrom' or 'meter':   #The program quits if it is passed inappropriate units (cf l 43, 55)
#         sys.exit('Please use Bohr, Angstrom, or meter. The program will now exit.')
      if self.units == 'Angstrom':
         self.geom /= 0.529177208 
         self.units = 'Bohr'
      elif self.units == 'Meter':
         self.geom /= ( 0.529177208 * ( 1 ** -10 ) )
         self.units = 'Bohr'
      else:
         print('Units are already in Bohrs.')

   def to_angstrom(self):
      #if self.units != 'Bohr' or 'Angstrom' or 'meter':
      #   sys.exit('Please use Bohr, Angstrom, or meter. The program will now exit.')
      if self.units == 'Bohr':
         self.geom *= 0.529177208
         self.units = 'Angstrom'
      elif self.units == 'Meter':
         self.geom *= ( 0.529177208 * (1 ** 10))
         self.units = 'Angstrom'
      else:
         print('Units are already in Angstroms.')

   def to_meter(self):
      #if self.units != 'Bohr' or 'Angstrom' or 'meter':
      #   sys.exit('Please use Bohr, Angstrom, or meter. The program will now exit.')
      if self.units == 'Angstrom':
         self.geom *= (1 ** -10)
         self.units = 'meter'
      elif self.units == 'Bohr':
         self.geom *= ( 0.529177208 * (1 ** -10))
         self.units = 'meter'
      else:
         print('Units are already in meters.')

   def xyz_string(self):
      string = '%d \n%s \n' %(self.natom, self.units)
      for i in range(self.natom):
         string += '%s%20.10f%20.10f%20.10f\n' %(self.labels[i],self.geom[i,0],self.geom[i,1],self.geom[i,2])
      return string
   
   def copy(self):
      copy = Molecule(self.xyz_string(),self.units)
      return copy

   def test_print(self):
      print(self.xyz)
      print(self.labels)
      print(self.masses)
      print(self.charges)
      print(self.geom)
      print(self.xyz_string())

mol = Molecule('molecule.xyz', 'Angstrom')
