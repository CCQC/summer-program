#! /usr/bin/python3

import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/aewiens')
from masses import get_mass
from molecule import Molecule
import numpy as np
from scipy import linalg as la


###################################### read in the moleulce
f = open("../extra-files/molecule.xyz").readlines()
mol = Molecule(f,"Bohr")
mol.angs()
atoms = mol.atoms
geom  = mol.geom
N = mol.__len__()


##################################### read in the Hessian matrix
g = open("../extra-files/hessian.dat","r")
H0 = np.matrix([i.split() for i in g.readlines()],float)


#################################### read in the atomic weights (amu) 
m = []                             # as a vector 
for i in atoms:
	m += [1/(get_mass(i))**0.5]*3
M = np.diag(m) 					     # as a matrix


################### mass-weight the Hessian & get eigenvalues/ vectors
mH = M*H0*M
e, l = la.eigh(mH)

########################  unit conversions
bohr2m = 5.2917721e-11
amu2kg = 1.6605389e-27
hartree2J = 4.3597443e-18
c = 29979245800.0                                              # cm /s
conv =  np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)


####################################   get frequencies and print to frequencies.txt
freq = [conv*np.sqrt(i) for i in e]
r = open("frequencies.txt","w")
for i in range(3*N):
	r.write("%12.7f\n" % freq[i])

###################################    visualize vibrational modes
t = open("modes.xyz","w")
for i in range(3*N):
	t.write("%d\n %s cm^{-1}\n" % (N, str(freq[i])))
	for j in range(N):
		t.write("%s%15.7f%15.7f%15.7f\n" % (atoms[j], geom[j][0], geom[j][1], geom[j][2]))
	t.write("\n")

l *= 0.52917721
print l
