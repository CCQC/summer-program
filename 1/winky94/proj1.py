#!/usr/bin/env python2
import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/winky94')
from molecule import Molecule
from masses import get_mass
import numpy as np

# Build molecule

mol = Molecule(open('../extra-files/molecule.xyz').read(), 'Bohr')
mol.bohr_to_ang()

# Read in hessian

lines = open('../extra-files/hessian.dat').readlines()

# Build mass-weighted hessian

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

# Compute eigenvalues and eigenvectors

e, v = np.linalg.eigh(mwH)

# un-mass-weight to normal coordinates

Q = np.array(np.matrix(M) * np.matrix(v))

# Get wavenumber

hartree2J = 4.3597443e-18
amu2kg = 1.6605389e-27
bohr2m = 5.2917721e-11
c = 29979245800.0 # speed of light in cm/s
convert = np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)

line_form = '{:2s}' + '{: >15.10f}'*6 + '\n'

out = ""
for a, eigh in enumerate(e):
    if eigh > 0:
        out += "{}\n{: >7.2f} cm^-1\n".format(len(mol), np.sqrt(eigh)*convert)
    else:
        out += "{}\n{: >7.2f}i cm^-1\n".format(len(mol), np.sqrt(-eigh)*convert)
    for i in range(len(mol)): 
        atom = mol.atoms[i]
        x, y, z = map(float,mol.coords[i])
        dx, dy, dz = Q[3*i:3*i + 3, a]
        out += line_form.format(atom, x, y, z, dx, dy, dz)
    out += '\n'

print(out)
