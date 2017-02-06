#!/usr/bin/python

import sys, numpy as np, conv as cnv
from molecule import Molecule

sys.path.insert(0,'../..')

class Frequencies:

	def __init__(self,mol,hessianString,G):

		self.mol  = mol
		self.hess = hessianString
		self.N    = mol.__len__()
		self.G    = G

		self.readHessian()

	def readHessian(self):

		#F  = np.matrix( [i.split() for i in self.hess.splitlines()],float )

		self.GF     = np.dot(self.G,F)
		self.Nmodes = self.GF.shape[0]
		
	def getFrequencies(self):

		self.e, self.l = np.linalg.eigh(self.GF)
		freq = []
		conv =  np.sqrt(cnv.hartree2J/(cnv.amu2kg*cnv.bohr2m**2))/(cnv.c*2*np.pi)  # dimensional analysis
		for i in self.e:
			if i <0:
				freq.append((-i)**0.5*conv)
			else:
				freq.append(i**0.5*conv)

		return freq

if __name__ == "__main__":

	f    = open("template.dat","r").read()
	mol  = Molecule(f)
	h    = open("hessian.dat","r").read()
	freq = Frequencies(mol,h,G)

	print( freq.getFrequencies() )
