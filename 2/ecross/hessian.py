#!/anaconda/bin/python

import sys
import numpy as np
from numpy import linalg as la 
import os

#imports the file for the Molecule clas
sys.path.insert(0, '../../0/ecross')
import molecule as M

#initiates a molecule from designated file with indicated units
mol = M.Molecule('./molecule.xyz', 'Bohr')

#file address for template to be used
templatefile =  './template.dat'

#creates template class to import molecule
class Template(object):
   def __init__(self,template,mol,command='psi4',disp_size=0.005,directory='DISPS'):
      """
      ===============================================================
      Construct a Hessian matrix from a set of molecular coordinates.
      ===============================================================
      
      Inputs:
      1) 'template.dat' file formatted consistent to the demands of 
         the chosen computational suite
      2) 'molecule.xyz' file containing molecular Cartesian coodinates
      3) the command to execute the computational software (string, default = 'psi4')
      4) size of displacements (float, default = 0.005)
      5) name of directory to contain displacements (string, default = 'DISPS')

      Returns:
      1) Directories, input.dat files, and output.dat files corresponding 
         to each of the displacements necessary to compute the Hessian 
         matrix
      2) Numpy Hessian matrix saved self.hessian and written to hessian.dat
      """
      self.mol = mol
      self.disp_size = disp_size
      self.directory = directory
      self.command = command

      File = open(template,'r')
      self.template = File.read() 
      File.close()

      self.generate_inputs()
      self.run_jobs()
      self.hessian = self.build_hessian()
      

   #saves mol parameters as 'input.dat'
   def write_file(self,coords):
      """
      Input: A template.dat file and a set of Cartesian coordinates for
      each atom in the molecule
      
      Output: Appropriately formatted file 'input.dat' containing units
      and the coordinates for each atom
      """
      xyz_string = ''
      for m in range(mol.natom):
         xyz_string += '{:s}{:17.12f}{:17.12f}{:17.12f}\n'.format(mol.labels[m],coords[3*m],coords[(3*m)+1],coords[(3*m)+2])
      a = self.template.format('units ' + str(mol.units) + '\n' + xyz_string)
      inputfile = open('input.dat','w')
      inputfile.write(a)
      inputfile.close()

   def generate_inputs(self):
      """
      Inputs: displacement size (float, default = 0.005)
              directory name (string, default = 'DISPS)

      Output: organized directories each containing 'input.dat' file
              according to the format in the indicated template
      """
      #create displacement directory, set as home_dir
      if not os.path.isdir(self.directory):
         os.mkdir(self.directory)
      home_dir = os.path.abspath(self.directory)
      os.chdir(home_dir)

      disps = [self.disp_size, -self.disp_size]
      
      """ NAUGHTH DISPLACEMENT """ 
      if not os.path.isdir('0_00'):
         os.mkdir('0_00')
      os.chdir('0_00')

      ref_coords = [ mol.geom[k,l] for k in range(mol.natom) for l in range(3) ]
      self.write_file(ref_coords)

      #return to displacement directory
      os.chdir(home_dir)
     
      """ FIRST DISPLACEMENTS """
      for i in range(mol.natom*3):
         for j in range(len(disps)):
            coords = [ mol.geom[k,l] for k in range(mol.natom) for l in range(3) ]
            coords[i] = coords[i] + disps[j]

            #make, enter directory for the displacement j in i
            b = '1_{:d}{:d}'.format(i,j)
            if not os.path.isdir(b):
               os.mkdir(b)
            os.chdir(b)
            
            #call the command to write input.dat
            self.write_file(coords)

            #return to displacement directory 
            os.chdir(home_dir)

      """ SECOND DISPLACEMENTS """ 
      for i in range(mol.natom*3):
         for j in range(i):
            for m in range(len(disps)):
               coords = [ mol.geom[k,l] for k in range(mol.natom) for l in range(3) ]
               coords[i] = coords[i] + disps[m]
               coords[j] = coords[j] + disps[m]

               b = '2_{:d}{:d}{:d}'.format(i,j,m)
               if not os.path.isdir(b):
                  os.mkdir(b)
               os.chdir(b)

               self.write_file(coords)

               #return to displacement directory
               os.chdir(home_dir)
      
      os.chdir('..')

   def run_jobs(self):
      """
      Inputs: self

      Output: 'output.dat' file from computational program for 
              each geometry in the displacement directory
      """
      d = os.listdir(self.directory)

      os.chdir(self.directory)

      for i in range(len(d)):
         os.chdir(d[i])
         os.system(self.command)
         os.chdir('..')

      os.chdir('..')
      

   def get_energy(self,match):
      """
      Search through Psi4 output.dat file and return the density-fitted
      Restricted Hartree-Fock final energy.
      """
      #sets an initial dummy value for energy
      energy = 0

      #opens output file 
      output = open('output.dat','r').read().splitlines()

      #search for energy from bottom of output; stop when energy is found
      for line in reversed(output):
         if line[:23] == match:
            energy = float(line.split()[-1])
            return(energy)
            break

   def build_hessian(self, energy_prefix = "  @DF-RHF Final Energy:"):
      os.chdir(self.directory)
      hessian = np.zeros((mol.natom*3,mol.natom*3))

      for i in range(mol.natom*3):
         for j in range(mol.natom*3):
            #computes diagonal elements of the Hessian
            if i == j:
               value = 0

               #subtracts two equivalents of the unperturbed energy from value
               os.chdir('0_00')
               value += - 2 * self.get_energy(energy_prefix)
               os.chdir('..')

               #adds E(Xa + h) and E(Xa - h) to value
               for k in range(2):
                  os.chdir('1_{:d}{:d}'.format(j,k))
                  energy = self.get_energy(energy_prefix)
                  value += energy
                  os.chdir('..')
                  
               #Divides value by the square of disp_size
               value = value * (self.disp_size ** -2)

               #sets diagonal elements of hessian to value
               hessian[i,j] = value

            #computes off-diagonal elements of the Hessian
            elif i > j:
               value = 0 

               #adds two equivalents of the unperturbed energy to value
               os.chdir('0_00')
               value += 2 * self.get_energy(energy_prefix)
               os.chdir('..')

               #subtracts E(Xa +/- h,Xb) and E(Xa, Xb +/- h) from value
               for k in range(2):
                  os.chdir('1_{:d}{:d}'.format(i,k))
                  energy = self.get_energy(energy_prefix)
                  value += -energy
                  os.chdir('..')

                  os.chdir('1_{:d}{:d}'.format(j,k))
                  energy = self.get_energy(energy_prefix)
                  value += - energy
                  os.chdir('..')
               
               #adds E(Xa +/- h, Xb +/- h) terms to value
               for l in range(2):
                  os.chdir('2_{:d}{:d}{:d}'.format(i,j,l))
                  energy = self.get_energy(energy_prefix)
                  value += energy
                  os.chdir('..')
              
               value = value * 0.5 * (self.disp_size ** -2)

               hessian[i,j] = value
               hessian[j,i] = value

      os.chdir('..')       
      output = np.savetxt('hessian.dat', hessian, '%20.15f',' ')            

      return(hessian)

if __name__ == '__main__':
   z = Template(templatefile,mol)
   print(z.hessian)
