import sys
from pprint import pprint
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/adabbott')
import masses
from molecule import Molecule
import numpy as np


mol = Molecule(open("../extra-files/molecule.xyz").read())

mol.bohr_to_ang()
#read Hessian Values

hessian_values = open('../extra-files/hessian.dat').readlines()
#collect hessain values into matrix

H = []
A = len(mol)  #this is redefined to give number of atoms
for row in hessian_values:
    for x in row.split():
        H.append(x)
hessian = np.matrix(H, float)
hessian = np.reshape(hessian,(3*A, 3*A))


m = []
M = []
for atom in mol.atoms:
    m.append(masses.get_mass(atom))
    m.append(masses.get_mass(atom))
    m.append(masses.get_mass(atom))
for mass in m:
    M.append(mass**-0.5)
M = np.diag(M)
mass_weighted_hess = M*hessian*M

eigens = np.linalg.eigh(mass_weighted_hess)
