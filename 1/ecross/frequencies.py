"""
TO DO:
-----

Fix the units problem (import file in Bohr screws with things)
Clean up script

"""
#!/usr/bin/python

import sys
import numpy as np

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

from numpy import linalg as la

evalues,evectors = la.eigh(mwhessian)

"""
Here we un-mass-weight the eigenvectors to get normal coordinates. 
"""

nevectors = np.dot(mass_mat, evectors)

"""
Here we determine the spatial frequencies and force constants for the H2O molecule
"""


