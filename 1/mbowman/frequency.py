#__author__ = "mbowman"

import sys
sys.path.insert(1, '../../extra-files')

import masses as m

import molecule

import numpy as np

import math

import cmath

class Frequency(object):

	def __init__(self, mole_str, hess_str):
		"""
		
		:param mole_str: string used to create instance of a new molecule
		:param hess_str: string used to create Hessian matrix
		"""
		self.units = "hertz"
		self.mol = molecule.Molecule(mole_str)
		self.mol.to_bohr()
		self.parse_hessian(hess_str)
		self.mwhess = self.mass_weight_mat(self.hess, self.mol.masses)
		self.eigVal, self.eigVec = np.linalg.eig(self.mwhess)
		"""	for a in range(len(self.eigVal)):
			print np.dot(self.mwhess, self.eigVec[a])
			print self.eigVal[a]*self.eigVec[a] """
		self.normCoord = np.array([ [(self.eigVec[b][a] / math.sqrt(self.mol.masses[b / 3])) for b in range(len(self.eigVec[a])) ] for a in range(len(self.eigVec)) ] )
		
		self.eigVal *= 4.35974465054 * 10**(-18) #hatree to joule
		self.eigVal /= (1.6605389 * 10**(-27) ) #amu to k
		self.eigVal /= ((5.2917721067 * 10**(-11) )**2) #bohr to m
		self.spFreq = np.empty([len(self.eigVal)], dtype=complex)
		for a in range(len(self.eigVal)):
			self.spFreq[a] = cmath.sqrt(self.eigVal[a]) / (2.0 * math.pi)
		self.to_wavenumber()
		with open('project-1-answers.xyz',"w+") as f:
			f.write(self.xyz_string())
		
	def parse_hessian(self, hess_str):
		"""

		:param hess_str: string used to create Hessian matrix
		"""			
		"""lines = hess_str.strip().split("\n")
		tempArray = [i.split("   ")[0:] for i in lines[0:]]
		for a in range(len(tempArray)):
			for b in range(len(tempArray[a])):
				tempArray[a][b] = float(tempArray[a][b])		"""
		self.hess = np.array([[float(b) for b in a.strip().split('   ')] for a in hess_str.splitlines()])

	def mass_weight_mat(self, mat, mass):
		"""
		Returns a mass weighted matrix from a unweighted matrix and a list of masses
		:param mat: 3N x 3N np float array for N atoms, in the case of a Hessian matrix, this describes the mixed partial derivatives
		:param mass: N float list of masses for atoms
		"""
		tempMatrix = np.array([ [(mat[a][b] / ( math.sqrt (mass[ a / 3] * mass[b / 3] ) ) )  for b in range(len(mat[a])) ] for a in range(len(mat)) ])
		return tempMatrix 

	def to_wavenumber(self):
		if self.units =="hertz":
			self.spFreq /= 29979245800
			self.units = "wavenumber"

	def to_hertz(self):
		if self.units == "wavenumber":
			self.spFreq *= 29979245800
			self.units = "hertz"
		 
	def xyz_string(self):
		str = ""
		for n in range(len(self.normCoord)):
			str += '{0:d} \n'.format(self.mol.natom) 
			if  self.spFreq[n].real >= self.spFreq[n].imag:
				str += '{0:6.2f} cm^-1\n'.format(self.spFreq[n].real)
			else:
				str += '{0:6.2f}i cm^-1\n'.format(self.spFreq[n].imag)
			for l in range(self.mol.natom):	
				str += self.mol.labels[l] + " " 
				str += '{0:16.10f}'.format(self.mol.geom[l][0])
				str += '{0:16.10f}'.format(self.mol.geom[l][1])
				str += '{0:16.10f}'.format(self.mol.geom[l][2])
				str += '{0:16.10f}'.format(self.normCoord[n][3*l+0])
				str += '{0:16.10f}'.format(self.normCoord[n][3*l+1])
				str += '{0:16.10f}'.format(self.normCoord[n][3*l+2])
				str += "\n"
			str += "\n"
		return str		
		

if __name__ == "__main__":
	with open("../../extra-files/hessian.dat") as openHess:
		hess_str = openHess.read()
	with open("../../extra-files/molecule.xyz") as openMol:
		mole_str = openMol.read()
	freq = Frequency(mole_str, hess_str)
