#! /usr/bin/python3

import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/aewiens')

import masses
from conv import bohr2m, amu2kg, hartree2J, c
from molecule import Molecule

import numpy as np

#####  Read in the molecule
f = open("../extra-files/molecule.xyz").read()
mol = Molecule(f,"Bohr")
mol.angs()
N = mol.__len__()


#####  Read in the Hessian
g = open("../extra-files/hessian.dat","r")
H0 = np.matrix([i.split() for i in g.readlines()],float)


#####  Construct M^{-1/2} (diagonal matrix)
m = []
for i in mol.atoms:
	m += [1/(masses.get_mass(i))**0.5]*3
M = np.diag(m)


#####  Mass-weight the Hessian
mH = M*H0*M


#####  Diagonalize Hessian
e, l = np.linalg.eigh(mH)


#####  Eigenvectors --> vibrational modes (in non-mass-weighted coordinates)
Q = np.matrix(M)*np.matrix(l)

#####   get frequencies
#####   lambda_a = omega_a^2 / conv^2
freq = []
conv =  np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)  # dimensional analysis
for i in e:
    if i <0:
        freq.append((-i)**0.5*conv)
    else:
        freq.append(i**0.5*conv)

print( freq )


#####    visualize vibrational modes
t = open("modes.xyz","w")
for i in range(3*N):
   t.write("%d\n %s cm^{-1}\n" % (N, str(freq[i])))
   for j in range(N):
      atom = mol.atoms[j]
      x,y,z = mol.geom[j,0], mol.geom[j,1], mol.geom[j,2]
      dx,dy,dz = Q[3*j,i], Q[3*j+1,i], Q[3*j+2,i]
      t.write("%s%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n" % (atom, x, y, z,dx,dy,dz))
   t.write("\n")

