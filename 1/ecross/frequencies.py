#!/anaconda/bin/python

import sys
import numpy as np
pi = np.pi
from numpy import linalg as la
import cmath 

sys.path.insert(0, '../../0/ecross') #Addres for the Molecule class file
import molecule as M

"""
The Hessian .dat file is converted to a numpy matrix, called 'hessian.'

Select Hessian file by changing the address 
found in hess_file
"""

hess_file = '../extra-files/hessian.dat'

hess_open = open(hess_file, 'r')
hess_string = hess_open.read()    #This creates a string of the Hessian
hess_split = hess_string.splitlines()     #This splits the string into a list of the lines in the string
hess_lol = []      
for i in hess_split:
   hess_lol.append(i.split())    #This splits the list of lines into a list of lists for the Hessian 

hessian = np.matrix(hess_lol, dtype=float)    #This converts the list of lists into a numpy matrix

"""
A Molecule object is built from the molecule.xyz file 
"""

sys.path.insert(0, '../extra-files') #Address for the molecule.xyz file
mol = M.Molecule('molecule.xyz','Bohr')

"""
Here we create the mass-weighted Hessian, labeled 'mwhessian.'
"""

s = 3 * mol.natom
a = (s,s)

mass_mat = np.zeros(a, dtype=float)   #Generates 3*natom x 3*natom matrix
molwt = []

for i in range(mol.natom):   #Generates 3*natom vector of the form [Ma Ma Ma Mb Mb Mb ... Mn Mn Mn]
   for j in range(3):
      molwt.append(mol.masses[int(i)])

for i in range(s):   #Generates a diagonal 3*natom square matrix with molwt along the diagonal
   mass_mat[i,i] = molwt[i] ** -0.5

mwhessian = np.dot((np.dot(mass_mat, hessian)),mass_mat)

"""
Computes the eigenvalues and eigenvestors of the mass-weighted Hessian matrix
"""

evalues,evectors = la.eigh(mwhessian)

"""
Here we un-mass-weight the eigenvectors to get normal coordinates. 
"""

nevectors = np.dot(mass_mat, evectors)
nevectors *= 0.529177208

"""
Here we determine the spatial frequencies and force constants for the H2O molecule
"""

hartree2j = ( 4.3597438e-18 )
bohr2m = ( 5.29177208e-11 ) 
amu2kg = ( 1.66054e-27 )
c = ( 2.99792458e10 ) 

#Returns Hessian eigenvalues in [rad(2) * s(-2)]
evalues_si = [(val * hartree2j / bohr2m / bohr2m / amu2kg ) for val in evalues]

#Returns frequencies in Hz 
vfreq_Hz = [1 / ( 2 * pi ) * (np.sqrt(np.complex_(val))) for val in evalues_si]

#Returns frequencies in cm(-1)
vfreq = [ ( val / c ) for val in vfreq_Hz ]


"""
Finally, we write the frequencies to a file in the .xyz format
"""
mol.to_angstrom()

string = ''

for w in range(len(nevectors)):
   string += '{:d}\n'.format(mol.natom)
   string += '{:9.2f} cm-1\n'.format(vfreq[w])
   for atom in range(mol.natom):
      string += '{:s}\t'.format(mol.labels[atom])
      string += '{:>15.10f}{:>15.10f}{:>15.10f}\t'.format(mol.geom[atom,0],mol.geom[atom,1],mol.geom[atom,2])
      string += '{:>15.10f}{:>15.10f}{:>15.10f}\n'.format(nevectors[w,3*atom],nevectors[w,3*atom+1],nevectors[w,3*atom+2])
   string += '\n'


f = open('output.xyz','w')
f.write(string)
f.close()





