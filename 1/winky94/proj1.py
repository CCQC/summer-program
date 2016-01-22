import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/winky94')
from molecule import Molecule
from masses import get_mass
import numpy as np

#read in the molecule
mol = Molecule(open('../extra-files/molecule.xyz').read(), 'Bohr')
mol.bohr_to_ang()

#read in the Hessian
lines = open('../extra-files/hessian.dat').readlines()
H = []
for line in lines:
    H.append(list(map(float, line.split())))
H = np.matrix(H)

#generate M^{-1/2}
M = []
for atom in mol.atoms:
    M.append(1/np.sqrt(get_mass(atom)))
    M.append(1/np.sqrt(get_mass(atom)))
    M.append(1/np.sqrt(get_mass(atom)))
M = np.diag(M)

#mass-weighted hessian
mwH = M * H * M

#mwH eigenvalues and vectors
E, V = np.linalg.eigh(mwH)

#vibrational modes
Q = np.matrix(M)*np.matrix(V)

#unit conversions
bohr2m = 5.2917721e-11
amu2kg = 1.6605389e-27
hartree2J = 4.3597443e-18
c = 29979245800.0                                              # speed of light, cm /s
conv =  np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)  # dimensional analysis

#vibratonal frequencies
freq = []
for e in E:
    if e < 0:
        freq.append((-e)**0.5*conv)
    else:
        freq.append(e**0.5*conv)

#genearte .xyz





print(mol)
