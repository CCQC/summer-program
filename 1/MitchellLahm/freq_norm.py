from Molecule import Molecule
import numpy as np
from numpy import linalg as LA
import math
import cmath

with open('hessian.dat','r') as my_file:
	data = my_file.read()
	data = data.split('\n')
	hessian = []
	data = data[:-1]
	for i in range(len(data)):
		temp = data[i].split()
		temp = map(float,temp)
		hessian.append(temp)
hessian = np.array(hessian)

Mol = Molecule('molecule.xyz')
Mol.to_angstrom()
def frequencies(Mol,hessian):
	mwhessian = hessian
	h = len(hessian)
	n = h/3
	M = np.identity(len(hessian))
	for i in range(n):
		for k in range(3):
			M[3*i+k] = M[3*i+k] / math.sqrt(Mol.masses[i])
	mwhessian = np.dot(mwhessian,M)
	mwhessian = np.transpose(mwhessian)
	mwhessian = np.dot(mwhessian,M)
	mwhessian = np.transpose(mwhessian)
	
	k, wq = LA.eig(mwhessian)
	q = np.dot(M,wq)
	q = np.transpose(q)
	
	n1 = 4.35974417e-18/((5.2917721092e-11**2)*1.6605389e-27)
	n2 = (2 * math.pi)
	c = 2.99792458e10

	k = [complex(i * n1,0) for i in k]
	v = np.sqrt(k)
	v = [i / n2 for i in v]
	v = np.array([i / c for i in v])
	
	# Sorting the v and q arrays in ascending vibrational frequency states
	v.imag = -v.imag
	l = np.array(np.argsort(v))
	v.imag = -v.imag
	q = q[l]
	v = v[l]
	xyz_string = ''
	for i in range(h):
		xyz_string += str(n) + '\n'
		if v[i].real == 0:
			xyz_string += '%7.2f'%(v[i].imag)+'i cm^-1\n'
		else:
			xyz_string += '%7.2f'%(v[i].real)+'  cm^-1\n'
		for j in range(n):
			xyz_string += '%s   %13.10f  %13.10f  %13.10f  %13.10f  %13.10f  %13.10f\n'%(Mol.labels[j],Mol.geom[j][0],Mol.geom[j][1],Mol.geom[j][2],q[i][3*j],q[i][3*j+1],q[i][3*j+2])
	
	return xyz_string

string = frequencies(Mol,hessian)
f = open('project1_answer.xyz','w')
f.write(string)
f.close()
			
			
			