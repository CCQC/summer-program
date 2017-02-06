#!/usr/bin/env python3

import numpy as np
import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/jevandezande')
from masses import get_mass
from molecule import Molecule

# Create a molecule
mol = Molecule(open('../extra-files/molecule.xyz').read(), 'Bohr')
mol.to_angstrom()

# Read the Hessian
H = np.genfromtxt('../extra-files/hessian.dat')

# Mass weight Hessian
# Build W = M^{-1/2}
weights = []
for atom in mol.atoms:
    w = 1/np.sqrt(get_mass(atom))
    weights += [w, w, w]
W = np.diag(weights)

# \Tilde H = M^{-1/2} H M^{-1/2}
Ht = W @ H @ W

# Diagonalize the Hessian
# \Tilde H = L \Lambda L^T
k, L = np.linalg.eigh(Ht)

# \lambda_a = \omega^2
hartree2J = 4.3597443e-18
amu2kg = 1.6605389e-27
bohr2m = 5.2917721e-11
c = 29979245800.0 # speed of light in cm/s
convert = np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)

# \Delta X(Q_sa) = Q_A M^{-1/2} l_A
Q = W @ L

out = ''
line_form = '{:2s}' + '{: >15.10f}'*6 + '\n'
for a, k_a in enumerate(k):
    if k_a < 0:
        out += '{}\n{: >7.2f}i cm^-1\n'.format(len(mol), np.sqrt(-k_a) * convert)
    else:
        out += '{}\n{: >7.2f}  cm^-1\n'.format(len(mol), np.sqrt(k_a) * convert)
    for i in range(len(mol)):
        atom = mol.atoms[i]
        x, y, z = mol.geom[i]
        dx, dy, dz = Q[3*i: 3*i + 3, a]
        out += line_form.format(atom, x, y, z, dx, dy, dz)
    out += '\n'

print(out)

