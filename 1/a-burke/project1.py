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
with open('../../extra-files/hessian.dat') as h:
	h = h.read()
h = h.split()
for i in range(len(h)):
	h[i] = float(h[i])
h = np.array(h)

with open('../../extra-files/molecule.xyz') as geometry:
	geometry = geometry.read()
temp = []
geometry = geometry.split()	
temp.append(geometry)


#Read in molecule.xyz
with open('../../extra-files/molecule.xyz') as geom:
	geom = geom.read()

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
	for i in range(len(con)):
		s = str(con[i])
		x += "\n"
		x += str(y) + "\n"	
		if s[-1] == 'j':
			s = s.replace('j','i')
			x += s + '\n'
		elif s[0] == '(':
			for n in s:
				s = s.replace('(','')
				s = s.replace(')','')
				s = s.replace('+0j','')
			x += s + '\n'
		val,val1,val2 = 0,1,2
		a,b,c,d = 0,3,4,5
		for j in range(len(mol.labels)):
			x += str(mol.labels[j]) + "  " + str(temp[a][b]) +   "  " + str(temp[a][c]) + "  " + str(temp[a][d]) + "  " +  str(q[i][val]) + "  " + str(q[i][val1]) + "  " + str(q[i][val2]) + '\n'
			b += 4
			c += 4
			d += 4
			val += 3
			val1 += 3
			val2 += 3
	
	file = open('output.dat', 'w')
	file.write(str(x))
	file.close()

frequencies(mol, h)
