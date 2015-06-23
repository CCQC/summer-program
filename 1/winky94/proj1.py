import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/winky94')
from molecule import Molecule
from masses import get_mass
import numpy as np

mol = Molecule(open('../extra-files/molecule.xyz').read(), 'Bohr')
mol.bohr_to_ang()

lines = open('../extra-files/hessian.dat').read().split('/n')

H = []

for line in lines:
    H.append(list(map(float, line.split())))

H = np.matrix(H)

M = []

for atom in mol.atoms:
    M.append(1/np.sqrt(get_mass(atom)))
    M.append(1/np.sqrt(get_mass(atom)))
    M.append(1/np.sqrt(get_mass(atom)))

M = np.diag(M)

mwH = M * H * M

print(mwH)
