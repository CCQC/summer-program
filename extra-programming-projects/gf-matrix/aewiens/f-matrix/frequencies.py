#!/usr/bin/python

import sys, numpy as np, conv as cnv
from molecule import Molecule

from masses import mass

##  build Gmatrix for water  ## 
muO  = 1/mass[8]
muH  = 1/mass[1]
r    = 1.09520/0.529177
phi  = 109*np.pi/180

G = np.zeros((3,3))
G[0,0] = G[1,1] = muO + muH 
G[1,0] = G[0,1] = muO*np.cos(phi)
G[0,2] = G[2,0] = -muO*np.sin(phi)/r
G[2,1] = G[1,2] = -muO*np.sin(phi)/r
G[2,2] = 2*(muO + muH - muO*np.cos(phi))/r**2

class Frequencies:

	def __init__(self,mol,hessianString,G):

		self.mol  = mol
		self.hess = hessianString.strip()
		self.N    = mol.__len__()
		self.G    = G
		self.readHessian()


	def readHessian(self):

		self.F      = np.matrix( [i.split() for i in self.hess.splitlines()],float )
		self.GF     = np.dot(self.G,self.F)
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

		return freq[::-1]


if __name__ == "__main__":

	f    = open("template.dat","r").read()
	mol  = Molecule(f)
	#h    = open("hessian.dat","r").read()
	h     = open("examples/h2o/FCMINT","r").read()
	freq = Frequencies(mol,h,G)

	print( freq.getFrequencies() )
