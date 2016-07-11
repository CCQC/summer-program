#!/anaconda/bin/python

import sys
import numpy as np
pi = np.pi
from numpy import linalg as la
import cmath 

sys.path.insert(0, '../../0/ecross')
import molecule as M

hess_file = '../extra-files/hessian.dat'

sys.path.insert(0, '../extra-files')
mol = M.Molecule('molecule.xyz','Bohr')

"""
=================================================================
Calculate the frequencies and normal modes from a Hessian matrix.
=================================================================

Inputs:
1) .xyz file containing molecule geometry in cartesian coordinates
2) .dat file containing 3N x 3N hessian matrix for specified molecule

Outputs:
1) File "output.xyz" containing the normal modes and their frequencies
2) A variety 
"""

class Frequencies(object):
   def __init__(self, molecule, hessian):
      self.mol = molecule
      self.hess_mat = self.read_hessian(hessian)
      self.mass_mat = self.mass_matrix()
      self.evalues, self.evectors = self.diagonalize()
      self.nmodes = self.unweight()
      self.freq = self.frequencies()
      self.write_file()

   """Read .dat file containing hessian matrix into a numpy matrix."""
   def read_hessian(self,hess_file):
      hess_string = open(hess_file, 'r').read()
      hess_split = hess_string.splitlines()
      hess_list = [ i.split() for i in hess_split ]      
      hess_mat = np.matrix( hess_list, dtype = float )
      return(hess_mat)

   """Create a 3N x 3N diagonal matrix to weight (unweight) the Hessian (eigenvectors)."""
   def mass_matrix(self):
      mass_mat = np.zeros( (3*self.mol.natom,3*self.mol.natom), dtype = float )
      molwt = [ mol.masses[int(i)] for i in range(self.mol.natom) for j in range(3) ]
      for i in range(len(molwt)):
         mass_mat[i,i] = molwt[i] ** -0.5 
      return(mass_mat)

   """Mass-weight the Hessian matrixi."""
   def mass_weight(self):
      mass_mat = self.mass_matrix()
      mwhessian = np.dot((np.dot(mass_mat,self.hess_mat)),mass_mat)
      return(mwhessian)
      
   """Diagonalize the mass-weighted hessian."""
   def diagonalize(self):
      mwhessian = self.mass_weight()
      evalues,evectors = la.eigh(mwhessian) 
      #Numpy is inscrutable and saves eigenvectors as columns #lame (hence transpose)
      evectors = np.transpose(evectors)
      return(evalues,evectors)

   """Unweight evectors, return normal modes."""
   def unweight(self):
      nmodes = np.dot(self.evectors,self.mass_mat) 
      return(nmodes)
      
   """Calculate frequencies from Hessian eigenvalues."""
   def frequencies(self):
      hartree2j = (4.3597438e-18)
      bohr2m = (5.29177208e-11) 
      amu2kg = (1.66054e-27)
      c =  (2.99792458e10) 
      evalues_si = [(val*hartree2j/bohr2m/bohr2m/amu2kg) for val in self.evalues]
      vfreq_hz = [1/(2*pi)*np.sqrt(np.complex_(val)) for val in evalues_si]
      vfreq = [(val/c) for val in vfreq_hz]
      return(vfreq)

   """Write file 'output.xyz' containing normal modes and their frequencies"""
   def write_file(self):
      mol.to_angstrom()
      string = ''
      for w in range(len(self.nmodes)):
         string += '{:d}\n'.format(mol.natom)
         string += '{:9.2f} cm-1\n'.format(self.freq[w])
         for atom in range(mol.natom):
            string += '{:s}\t'.format(mol.labels[atom])
            string += '{:>15.10f}{:>15.10f}{:>15.10f}\t'.format(mol.geom[atom,0],mol.geom[atom,1],mol.geom[atom,2])
            string += '{:>15.10f}{:>15.10f}{:>15.10f}\n'.format(self.nmodes[w,3*atom],self.nmodes[w,3*atom+1],self.nmodes[w,3*atom+2])
         string += '\n'
      with open('output.xyz','w') as f: 
         f.write(string)

if __name__ == '__main__':
   test = Frequencies(mol, hess_file)
   
