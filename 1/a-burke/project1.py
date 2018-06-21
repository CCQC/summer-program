#Author: Ally Burke

import sys
sys.path.append('../../extra-files')
sys.path.append('../../0/a-burke')
from molecule import Molecule
import numpy as np
from numpy import linalg as la
from masses import mass
from masses import charge
import numpy.lib.scimath as nls
import scipy.constants as sci
import math

#Read in hessian matrix
h = open('../../extra-files/hessian.dat').read()
h = h.split()
for i in range(len(h)):
    h[i] = float(h[i])
h = np.array(h)

geometry = open('../../extra-files/molecule.xyz').read()
temp = []
geometry = geometry.split()	
temp.append(geometry)


#Read in molecule.xyz
geom = open('../../extra-files/molecule.xyz').read()

#Build a Molecule object from molecule.xyz
mol = Molecule(geom, "Angstrom")

def frequencies(mol, h):
	
	#reshapes the hessian to 9 3x3 matrices
	N = mol.natom
	h = h.reshape(3*N, 3*N)
		
	#array of masses
	m = []	
	
	#Build the mass-weighted hessian matrix	
	for y in mol.labels:
		m.append([mass[charge[y]]]*3) 
	m = np.sqrt(m)	
	M = np.diagflat(m)
	M = la.inv(M)
	H = np.dot(M,h) 
	H = np.dot(H,M)
	
	#Compute the eigenvalues and eigenvectors of the mass-weighted Hessian
	#K = eigenvalues, Q = eigenvectors
	K,Q = np.linalg.eigh(H)
	
	#Un-mass-weight the eigenvectors
	q = np.dot(H,Q)
	q = np.transpose(q)
	#Determine the spatial frequencies from K in wavenumbers (cm-1)
	#Conversion factors
	k = (K*sci.physical_constants['hartree-joule relationship'][0])/(((sci.physical_constants['atomic unit of length'][0])**2)*(sci.physical_constants['atomic mass unit-kilogram relationship'][0])) #converting to J, m2, and kg
	cn = nls.sqrt(k) / (2*math.pi)
	con = cn/(sci.physical_constants['speed of light in vacuum'][0]*100)
			
	
	#Write out the normal modes to output.dat in .xyz format
	#x is written to the output file
	x=''
	q = np.array(q)
	y = mol.natom
	b = 3
	c = 4
	d = 5
	num = 0
	val = 0
	val1 = 1
	val2 = 2
	for i in range(len(con)):
		x += "\n"
		x += str(y) + "\n"	
		x += '{:f}\n'.format(con[i])
		val = 0
		val1 = 1
		val2 = 2
		a = 0
		b = 3
		c = 4
		d = 5
		for j in range(len(mol.labels)):
			x += str(mol.labels[j]) + "  "
			x += str(temp[a][b]) +   "  " + str(temp[a][c]) + "  " + str(temp[a][d]) + "  " +  str(q[num][val]) + "  " + str(q[num][val1]) + "  " + str(q[num][val2])
			b += 4
			c += 4
			d += 4
			val += 3
			val1 += 3
			val2 += 3
			x += '\n'
		num += 1
	
	file = open('output.dat', 'w')
	file.write(str(x))
	file.close()

frequencies(mol, h)
