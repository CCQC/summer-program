#!/usr/bin/env python2

import sys
sys.path.insert(0, '../extra-files')
from masses import get_mass
import numpy as np
import scipy.linalg as la
from molecule import Molecule
from pprint import pprint

# Read molecule
mol = Molecule(open('../extra-files/molecule.xyz').read(), 'Bohr')
mol.to_angstrom()

# Read the Hessian
lines = open('../extra-files/hessian.dat').readlines()

H = []
for line in lines:
    H.append(list(map(float, line.split())))

H = np.matrix(H)

# Mass weight Hessian
# \Tilde H = M^{-1/2} H M^{-1/2}

# Build M^{-1/2}
masses = []
for atom in mol.atoms:
    masses.append(1/np.sqrt(get_mass(atom)))
    masses.append(1/np.sqrt(get_mass(atom)))
    masses.append(1/np.sqrt(get_mass(atom)))
M = np.diag(masses)

Ht = M * H * M

# Diagonalize the Hessian
# \Tilde H = L \Lambda L^T

k, L = la.eigh(Ht)

# \lambda_a = \omega^2

hartree2J = 4.3597443e-18
amu2kg = 1.6605389e-27
bohr2m = 5.2917721e-11
c = 29979245800.0 # speed of light in cm/s
convert = np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)

# \Delta X(Q_sa) = Q_A M^{-1/2} l_A
Q = np.matrix(M) * np.matrix(L)

out = ''
line_form = '{:2s}' + '{: >15.10f}'*6 + '\n'
for a, k_a in enumerate(k):
    if k_a < 0:
        out += '{}\n{: >7.2f}i cm^-1\n'.format(mol.num, np.sqrt(-k_a) * convert)
    else:
        out += '{}\n{: >7.2f}  cm^-1\n'.format(mol.num, np.sqrt(k_a) * convert)
    for i in range(mol.num):
        atom = mol.atoms[i]
        x, y, z = mol.geom[i]
        dx = Q[3*i, a]
        dy = Q[3*i + 1, a]
        dz = Q[3*i + 2, a]
        out += line_form.format(atom, x, y, z, dx, dy, dz)
    out += '\n'

print(out)

