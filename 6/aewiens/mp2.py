#!/usr/bin/env python3

import psi4, sys, numpy as np, configparser as cfp
sys.path.insert(0,"../../5/aewiens/")
from uhf import UHF

class MP2:

	def __init__(self,mol,mints,maxiter,conv):

		uhf = UHF(mol,mints,maxiter,conv)
		uhf.computeEnergy()

		self.nocc = uhf.nelec
		self.norb = uhf.norb
		self.E0   = uhf.E
		self.e    = uhf.e
		self.C    = uhf.C
		self.G    = uhf.G


	def tei_einsum(self,g,C):
		return np.einsum("Pp,Pqrs->pqrs",C,
				np.einsum("Qq,PQrs->Pqrs",C,
				np.einsum("Rr,PQRs->PQrs",C,
				np.einsum("Ss,PQRS->PQRs",C,g))))


	def computeEnergy(self):
		"""
		Spin-orbital implementation of mp2 equations
		"""

		Gmo  = self.tei_einsum(self.G, self.C)
		nocc = self.nocc
		norb = self.norb

		e  = self.e
		Ec = 0.0
		for i in range(nocc):
			for j in range(nocc):
				for a in range(nocc,norb):
					for b in range(nocc,norb):
						#Ecorr += (2*Gmo[i,j,a,b]-Gmo[i,j,b,a])*Gmo[i,j,a,b]/(e[i]+e[j]-e[a]-e[b])
						Ec += (Gmo[i,j,a,b]*Gmo[a,b,i,j])/(e[i]+e[j]-e[a]-e[b])

		Ec *= 0.25
						
		"""
		o = slice(0,nocc)
		v = slice(nocc,nocc+norb)
		x = np.newaxis

		D = e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v]
		T = Gmo[o,o,v,v]*Gmo[v,v,o,o] 
		#T =  np.square(Gmo[o,o,v,v])
		#T /= D

		Ecorr = 0.25*np.ndarray.sum(T,axis=(0,1,2,3))
		"""

		return self.E0 + Ec



if __name__ == '__main__':
	
	config = cfp.ConfigParser()
	config.read('Options.ini')

	molecule = psi4.geometry( config['DEFAULT']['molecule'] )
	molecule.update_geometry()

	basis = psi4.core.BasisSet.build(molecule, "BASIS", config['DEFAULT']['basis'],puream=0)
	mints = psi4.core.MintsHelper(basis)

	maxiter   = int( config['SCF']['maxIter'] )
	conv      = float( config['SCF']['conv']  )

	mp2 = MP2(molecule,mints,maxiter,conv)
	print( mp2.computeEnergy() )
