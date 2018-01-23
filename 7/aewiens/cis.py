#!/usr/bin/env python3
import sys, configparser, numpy as np
from scipy import linalg as la
sys.path.insert(0,"../../5/aewiens")
sys.path.insert(0,"../../6/aewiens")
from uhf import UHF
from mp2 import MP2

class CIS:

	def __init__(self,options):

		uhf = UHF( options )

		uhf.computeEnergy()

		self.nocc  = uhf.nocc
		self.nvirt = uhf.norb - self.nocc
		self.E0    = uhf.E
		self.e     = uhf.e

		mp2 = MP2( options )
		self.GMO   = mp2.transformTEI( uhf.G, uhf.C )

		print( "-----------------------------------" )
		print( "| Output for CI Singles Procedure |" )
		print( "-----------------------------------" )
		print( "\n @UHF Energy:   %f" % uhf.E )


	def getSingles(self):
		"""  return a list of all (i,a) single excitations  xi_i -> xi_a
		"""
		return [(i,a) for i in range(self.nocc) for a in range(self.nocc,self.nocc+self.nvirt)]


	def computeStates(self):

		E0 = self.E0
		e  = self.e
		nDeterminants = self.nocc*self.nvirt
		excitations   = self.getSingles()

		##  build CIS hamiltonian (tH)  ##
		tH = np.zeros((nDeterminants,nDeterminants))
		for P, (i,a) in enumerate(excitations):
			for Q, (j,b) in enumerate(excitations):
				tH[P,Q] = self.GMO[a,j,i,b] + (e[a] - e[i])*(a==b)*(i==j) 

		E, C = np.linalg.eigh(tH)

		self.E = E
		eigenvectors = C.T

		info = []
		for i, row in enumerate( C.T ):
			temp =  ""
			for j, value in enumerate(row):
				if value**2 >= 0.10:
					percentage = "{:5.0f}% ".format( 100*value**2 )
					excitation = " {:d} --> {:d}".format(*excitations[j])
					temp += percentage + excitation + "   "

			info.append( temp )

		print("\n State      Energy (Eh)    Excitations")
		print("-----------------------------------------------------------")
		for i, energy in enumerate(E):
			print("{:4d}  {: >16.11f}  ".format(i,energy) + info[i] )
		print("-----------------------------------------------------------")



if __name__ == '__main__':
	
	config = configparser.ConfigParser()
	config.read('Options.ini')

	cis = CIS(config)
	cis.computeStates()

