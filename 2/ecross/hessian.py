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

#creates template class to import molecule
class Template(object):
   def __init__(self,template,mol):
      self.mol = mol

      File = open(template,'r')
      self.template = File.read() 
      File.close()

   #saves mol parameters as 'input.dat'
   def write_file(self,coords):
      xyz_string = ''
      for m in range(mol.natom):
         xyz_string += '{:s}{:17.12f}{:17.12f}{:17.12f}\n'.format(mol.labels[m],coords[3*m],coords[(3*m)+1],coords[(3*m)+2])
      a = self.template.format('units ' + str(mol.units) + '\n' + xyz_string)
      inputfile = open('input.dat','w')
      inputfile.write(a)
      inputfile.close()

   def generate_inputs(self,disp_size=0.005,directory='DISPS'):
      
      #creates DISPS directory
      if not os.path.isdir(directory):
         os.mkdir(directory)
      os.chdir(directory)

      disps = [disp_size, -disp_size]
      #stores the reference geometry
      ref_coords = [ mol.geom[k,l] for k in range(mol.natom) for l in range(3) ]
      
      os.mkdir('0_00')
      os.chdir('0_00')

      self.write_file(ref_coords)
      os.chdir('..')
     
      #computes the first displacement for all coordinates of all atoms in both directions
      for i in range(mol.natom*3):
         for j in range(len(disps)):
            coords = [ mol.geom[k,l] for k in range(mol.natom) for l in range(3) ]
            coords[i] = coords[i] + disps[j]

            #makes, enters directory for the displacement j in i
            b = '1_{:d}{:d}'.format(i,j)
            os.mkdir(b)
            os.chdir(b)
            
            #calls the command to write input.dat
            self.write_file(coords)

            #returns to DISPS directory 
            os.chdir('..')

      for i in range(mol.natom*3):
         for j in range(i):
            for m in range(len(disps)):
               coords = [ mol.geom[k,l] for k in range(mol.natom) for l in range(3) ]
               coords[i] = coords[i] + disps[m]
               coords[j] = coords[j] + disps[m]

               b = '2_{:d}{:d}{:d}'.format(i,j,m)
               os.mkdir(b)
               os.chdir(b)

               self.write_file(coords)

               os.chdir('..')
      
      os.chdir('..')

   def run_jobs(self, command = 'psi4', directory = 'DISPS'):
      d = os.listdir(directory)

      os.chdir(directory)

      for i in range(len(d)):
         os.chdir(d[i])
         os.system(command)
         os.chdir('..')

      os.chdir('..')

   def get_energy(self,match = "  @DF-RHF Final Energy:"):
      #sets an initial dummy value for energy
      energy = 0

      #opens psi4 output, 
      output = open('output.dat','r').read().splitlines()

      #searches for energy from bottom of output, stops when energy is found
      for line in reversed(output):
         if line[:23] == match:
            energy = float(line.split()[-1])
            return(energy)
            break

   def build_hessian(self, energy_prefix = "  @DF-RHF Final Energy:", disp_size=0.005, directory='DISPS'):
      os.chdir(directory)
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
               value = value * (disp_size ** -2)

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
                  value += - energy
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
              
               value = value * 0.5 * (disp_size ** -2)

               hessian[i,j] = value
               hessian[j,i] = value

      os.chdir('..')       
      output = np.savetxt('hessian.dat', hessian, '%20.15f',' ')            
      return(output)

#file address for template to be used
templatefile =  './template.dat'
z = Template(templatefile,mol)

"""
z.generate_inputs()
z.run_jobs()
"""
z.build_hessian()

