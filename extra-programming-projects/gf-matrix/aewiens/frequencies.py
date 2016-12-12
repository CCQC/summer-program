#!/usr/bin/python

import numpy as np
from molecule import Molecule

##   Conversions   ##

bohr2m    = 5.2917721e-11
amu2kg    = 1.6605389e-27
hartree2J = 4.3597443e-18
c         = 29979245800.0

class Frequencies:

	def __init__(self,mol,hessianString):

		self.mol = mol
		self.hessianString = hessianString
		self.N = mol.__len__()

		m = []
		for i in range(self.N):
			m += [1/(mol.masses[i])**0.5]*3
		self.MM = np.diag(m)
		self.m = m

	def get_MWhessian(self):

		H0 = np.matrix( [i.split() for i in self.hessianString.splitlines()],float )
		mwH = np.dot(self.MM, np.dot(H0,self.MM) )
		return mwH

	def get_frequencies(self):

		self.e, self.l = np.linalg.eigh( self.get_MWhessian() )
		self.Q = np.matrix(self.MM)*np.matrix(self.l)
		freq = []
		conv =  np.sqrt(cnv.hartree2J/(cnv.amu2kg*cnv.bohr2m**2))/(cnv.c*2*np.pi)  # dimensional analysis
		for i in self.e:
			if i <0:
				freq.append((-i)**0.5*conv)
			else:
				freq.append(i**0.5*conv)

		return freq

	def output(self, output):

		mol = self.mol
		freq = self.get_frequencies()

		t = open(output,"w")
		for i in range(3*self.N):
			t.write("%d\n%s cm^{-1}\n" % (self.N, str(freq[i])))
			for j in range(self.N):
				atom = mol.atoms[j]
				x,y,z = mol.geom[j,0], mol.geom[j,1], mol.geom[j,2]
				dx,dy,dz = self.Q[3*j,i], self.Q[3*j+1,i], self.Q[3*j+2,i]
				t.write("{:s}{:12.7f}{:12.7f}{:12.7f}\n".format(atom, x,y,z) )
			t.write("\n")

		return None
        
if __name__ == "__main__":
	f   = open("h2o.zmat","r").read()
	mol = Molecule(f,"Bohr")
	hessianString = open("hessian.dat","r").read()
	frequencyObj  = Frequencies(mol,hessianString)
	frequencyObj.output("modes.xyz")
