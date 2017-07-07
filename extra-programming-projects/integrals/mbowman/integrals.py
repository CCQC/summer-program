#__author__ = "mbowman"

import sys
sys.path.insert(0, '../../../extra-files')
sys.path.insert(0, '../../../0/mbowman')

import STO3Gbasis as sto
import molecule
import masses
import math
import numpy as np
from scipy import special, misc
import collections

class Integrals(object):
	
	def __init__(self, mol_str):
		"""
		""" 
		self.mol = molecule.Molecule(mol_str)
		self.mol.to_bohr()
		self.atoms = self.mol.labels
		self.coord = np.array(self.mol.geom)
		self.createBasisDict() 
		self.K = len(self.bd) #dimension of molecular integrals		
		self.gfn = 3 #number of Gaussian functions to sum over
		self.overlapIntegral()
		
	def subshells(self, z):
		s = ['1_s']
		if z > 2:
			s.extend(['2_s','2_p_x', '2_p_y', '2_p_z'])
			if z > 10:
				s.extend(['3_s','3_p_x', '3_p_y', '3_p_z']) 
		return s

	def createBasisDict(self):
		self.bd = {}
		duplicates = [item for item, count in collections.Counter(self.atoms).items() if count > 1]
		duplicateCounter = [1] * len(duplicates)
		for i, atom in enumerate(self.atoms):
			tempCharge = masses.get_charge(atom) 	#atomic number of atom, determines subshells
			if atom in duplicates:
				tas = atom + "_" + str(duplicateCounter[duplicates.index(atom)]) 	#creates a string to denote duplicate atoms e.g. H_2
				duplicateCounter[duplicates.index(atom)] += 1
			else:
				tas = atom							#creates a string to denote unique atoms e.g. O
			tempShells = self.subshells(tempCharge) 	#creates temp array of subshells from atom
			for subShell in tempShells:
				shellInfo = subShell.split('_') 	#splits shell into three sections based on three quantum numbers
				sss = tas + '_' + subShell
				self.bd[sss] = {} 		#initializes dictionary entry for subshell of atom
				self.bd[sss]['a'] = sto.expFact(tempCharge, int(shellInfo[0])) #assigns expansion factors based on atomic number and shell
				if shellInfo[1] == 's':
					self.bd[sss]['d'] = sto.contCoeff(int(shellInfo[0]),1)
					self.bd[sss]['l'] = 0
					self.bd[sss]['m'] = 0
					self.bd[sss]['n'] = 0
				elif shellInfo[1] == 'p':
					self.bd[sss]['d'] = sto.contCoeff(int(shellInfo[0]),2)
					if shellInfo[2] == 'x':
						self.bd[sss]['l'] = 1
						self.bd[sss]['m'] = 0
						self.bd[sss]['n'] = 0
					elif shellInfo[2] == 'y':
                                        	self.bd[sss]['l'] = 0
                                        	self.bd[sss]['m'] = 1
                                        	self.bd[sss]['n'] = 0
					elif shellInfo[2] == 'z':
						self.bd[sss]['l'] = 0
						self.bd[sss]['m'] = 0
						self.bd[sss]['n'] = 1 
					else:
						print "Error: p orbital not specified correctly"
				elif shellInfo[1] == 'd':
					print "Error: d orbitals not implemented"
				self.bd[sss]['R'] = [self.coord[i][c] for c in range(3)] 
	
	def overlapIntegral(self):			 
		self.S = np.zeros((self.K, self.K))
		for i, orb1 in enumerate(self.bd):
			for j, orb2 in enumerate(self.bd):
				for alphaA, da in zip(self.bd[orb1]['a'], self.bd[orb1]['d']):
					for alphaB, db in zip(self.bd[orb2]['a'], self.bd[orb2]['d']):
						gamma = alphaA + alphaB
						p = self.p(self.bd[orb1]['R'], self.bd[orb2]['R'], alphaA, alphaB)
						#print p
						t = da * db * self.normalization(alphaA, self.bd[orb1]['l'], self.bd[orb1]['m'],self.bd[orb1]['n'] ) * self.normalization(alphaB, self.bd[orb2]['l'], self.bd[orb2]['m'], self.bd[orb2]['n'])
						#print t
						norm = ( self.bd[orb1]['R'][0] - self.bd[orb2]['R'][0] + self.bd[orb1]['R'][1] - self.bd[orb2]['R'][1] + self.bd[orb1]['R'][2] - self.bd[orb2]['R'][2])
						#print norm
						t *= math.exp( -1*alphaA*alphaB*norm / gamma  )
						print t
						t *= self.f2(0 , gamma, self.bd[orb1]['l'], self.bd[orb2]['l'], self.bd[orb1]['R'], self.bd[orb2]['R'], p) 
						print t
						t *= self.f2(1 , gamma, self.bd[orb1]['m'], self.bd[orb2]['m'], self.bd[orb1]['R'], self.bd[orb2]['R'], p)
						print t
						t *= self.f2(2 , gamma, self.bd[orb1]['n'], self.bd[orb2]['n'], self.bd[orb1]['R'], self.bd[orb2]['R'], p)	
						print t
						print " "
		#return S	

	def f1(self, j, l, m, a, b):
		"""
		corresponds to equation 8
		"""
		s = 0
		for k in range(max(0,(j-m)), (min(j, l) + 1)):
			s += special.binom(l,k)*special.binom(m,(j-k))*(a ** (l - k) )*(b ** (m + k - j) )
		return s

	def p(self, a, b, alphaA, alphaB):
		"""
		corresponds to equation 5
		"""
		return ( np.multiply(alphaA, a) + np.multiply(alphaB, b) )/ (alphaA * alphaB)

	def f2(self, dimension, gamma, l, m, a, b, p):
		"""
		corresponds to equation 3 in notes
		""" 
		s = 0
		for x in range(int( (l + m)/2 )):
			s += self.f1(x, l, m, np.subtract(p,a)[dimension], np.subtract(p,b)[dimension]) * misc.factorial2((2*x-1),exact=True) / ((2* gamma) ** x) 
		s *= math.sqrt(( math.pi / gamma) )
		return s	
	
	def normalization(self, alpha, l, m, n):
		"""
		corresponds to equation 9
		"""
		s = (4 * alpha) ** (l + m + n)
		s /= ( misc.factorial2( (2*l -1), exact=True) * misc.factorial2( (2*m -1), exact=True) * misc.factorial2( (2*n -1), exact=True) )
		s *= math.pow( (2*alpha/math.pi), 1.5) 
		return math.sqrt(s)			
if __name__ == "__main__":
	mole_str = open("../../../extra-files/molecule.xyz").read()
	integral = Integrals(mole_str)
	print integral.S
