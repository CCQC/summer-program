#!/anaconda/bin/python 

import sys 
import numpy as np

class Molecule(object):
   def __init__(self, molecule, units):
      self.molecule = molecule
      self.units = units.lower()
      self.natom, self.labels, self.geom = self.read()
      self.masses, self.charges = self.mass_charge()
   
   def __len__(self):
      return self.natom

   """Read input file and return number of atoms, atom labels, and geometry."""
   def read(self):
      mol_string = open(self.molecule, 'r').read()
      mol_list = mol_string.splitlines()
      natom = int(mol_list[0])
      xyz = [i.split() for i in mol_list[2:natom+2] ]
      labels = [ i[0] for i in xyz ]
      geom = np.matrix( [ [ float(i[j]) for j in range(1,4) ] for i in xyz ] )
      return( natom, labels, geom ) 
  
   """Import masses/charges program and return masses and charges for molecule."""
   def mass_charge(self):
      sys.path.insert(0,'../../extra-files')
      import masses as M

      mass = [ float(M.get_mass(self.labels[i])) for i in range(self.natom) ]
      charge = [ M.get_charge(self.labels[i]) for i in range(self.natom) ]
      return( mass, charge )

   """
   Unit Conversion Functions:
   --------------------------
   Convert units to_bohr, to_angstrom, or to_meter.
   
   Arguments: self.units (must be bohr, angstrom, or meter)
   Change coordinate units in self.geom and string value of self.units.
   """
   def to_bohr(self):   
      if self.units == 'bohr':
         print('Units are already in Bohr.')
      elif self.units == 'angstrom':
         self.geom /= 0.529177208 
         self.units = 'bohr'
      elif self.units == 'meter':
         self.geom /= ( 0.529177208 * ( 1 ** -10 ) )
         self.units = 'bohr'
      else: 
         print('Warning: Molecule built with unacceptable units; conversion not processed.')

   def to_angstrom(self):
      if self.units == 'angstrom':
         print('Units are already in Angstroms.')
      elif self.units == 'bohr':
         self.geom *= 0.529177208
         self.units = 'angstrom'
      elif self.units == 'meter':
         self.geom *= ( 0.529177208 * (1 ** 10))
         self.units = 'angstrom'
      else:
         print('Warning: Molecule built with unacceptable units; conversion not processed.')
      
   def to_meter(self):
      if self.units == 'meter':
         print('Units are already in meters.')
      elif self.units == 'angstrom':
         self.geom *= (1 ** -10)
         self.units = 'meter'
      elif self.units == 'bohr':
         self.geom *= ( 0.529177208 * (1 ** -10))
         self.units = 'meter'
      else:
         print('Warning: Molecule built with unacceptable units; conversion not processed.')
   
   """Return atom labels and coordinates for molecule."""
   def geom_string(self):
      string = ''
      for i in range(self.natom):
         string += '{:s}{:20.10f}{:20.10f}{:20.10f}\n'.format(self.labels[i],self.geom[i,0],self.geom[i,1],self.geom[i,2])
      return string[:-1]

   """Return Molecule object as string in .xyz format."""
   def xyz_string(self):
      string = '{:d} \n{:s} \n'.format(self.natom, self.units)
      for i in range(self.natom):
         string += '{:s}{:20.10f}{:20.10f}{:20.10f}\n'.format(self.labels[i],self.geom[i,0],self.geom[i,1],self.geom[i,2])
      return string
   
   """Return a fresh copy of self."""
   def copy(self):
      inp = open('copy.xyz','w')
      inp.write(self.xyz_string())
      inp.close()
      copy = Molecule('copy.xyz',self.units)
      return copy


   
if __name__ == '__main__':
   mol = Molecule('molecule.xyz', 'Bohr')

   print('Atoms: ' + str(mol.labels) + '\n')
   print('Masses: ' + str(mol.masses) + '\n')
   print('Charges: ' + str(mol.charges) + '\n')
   print('Coordinates:\n' + str(mol.geom) + '\n')
   print('Output:\n\n' + mol.xyz_string())

