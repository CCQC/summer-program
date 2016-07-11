import sys
sys.path.insert(0, '../../extra-files')
sys.path.insert(0, '../../0/mjames84')
from molecule import Molecule
from masses import get_mass
import numpy as np

# Read hessian matrix
  
lines = open('../extra-files/hessian.dat').readlines()

# Build molecule object

geom_str = (open('../extra-files/molecule.xyz')).read()
mol = Molecule(geom_str)

# Build mass-weighted hessian 

hess = [] 
for line in lines:
    hess.append(list(map(float,line.split())))
hess = np.matrix(hess)


mass = []
for atom in mol.atoms:
    mass.append(1/np.sqrt(get_mass(atom)))
    mass.append(1/np.sqrt(get_mass(atom)))
    mass.append(1/np.sqrt(get_mass(atom)))
mass = np.diag(mass) 
mass_weighted_hess = mass * hess * mass

# Compute eigenvalues and eigenvectors of mass-weighted hess 

eigenval, eigenvect = np.linalg.eigh(mass_weighted_hess)

# Un-mass-weight eigenvectors to get normal coords

Q = np.array(np.matrix(mass) * np.matrix(eigenvect))

# Determine spatial frequencies (in cm-1) from force constant

#unit conversions

hartree2J = 4.3597443e-18
amu2kg = 1.6605389e-27
bohr2m = 5.2917721e-11
c = 29979245800.0 # cm/s
convert = np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)
freq = ""
line_form = '{:2s}' + '{: >15.10f}' + '\n'
for a, e in enumerate(eigenval):
    if e < 0:
        freq += "{}\n{: >7.2f}i cm^-1\n".format(len(mol), np.sqrt(-e)*convert)
    else:
        freq += "{}\n{: >7.2f}i cm^-1\n".format(len(mol), np.sqrt(e)*convert)
    for i in range(len(mol)):
        atom = mol.atoms[i]
        x, y, z = mol.xyz[i]
        dx, dy, dz = Q[3*i: 3*i + 3, a]
        freq += line_form.format(atom, x, y, z, dx, dy, dz)

print(freq)
 
