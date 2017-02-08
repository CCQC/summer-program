#!/usr/bin/env python3

import psi4, numpy as np, configparser as cfp

class UHF:

	def __init__(self,mol,mints,maxiter,conv):

		self.getIntegrals(mints)

		self.nelec   = self.getNelec(mol)
		self.conv    = conv
		self.maxiter = maxiter
		self.norb    = len(self.S)

		self.Vnu = mol.nuclear_repulsion_energy()
		self.E   = 0.0


	def getIntegrals(self,mints):

		self.T = block_oei( mints.ao_kinetic() )
		self.V = block_oei( mints.ao_potential() )
		self.S = block_oei( mints.ao_overlap() )                      

		G = block_tei(np.array( mints.ao_eri() ) )
		self.g = G.transpose((0,2,1,3))-G.transpose((0,2,3,1))

		S = mints.ao_overlap()
		S.power(-0.5,1.e-16)
		self.X = block_oei( S.to_array() )
		self.D = np.zeros_like(self.S)
		

	def getNelec(self, mol):

		char = mol.molecular_charge()
		nelec = -char
		for A in range(mol.natom()):
			nelec += mol.Z(A)

		return int(nelec)


	def computeEnergy(self):

		H = self.T + self.V
		X = self.X
		g = self.g
		D = self.D

		for i in range(self.maxiter):

			v = np.einsum("mnrs,ns->mr", g, self.D)
			F = H + v
			e,tC = np.linalg.eigh(X@F@ X)

			C  = X@tC
			oC = C[:,:self.nelec]
			D  = oC * oC.T

			E0 = self.E
			E  = np.trace( (H+0.5*v)*D) + self.Vnu
			dE = np.fabs(E-E0)
			
			if __name__ == '__main__':
				print("UHF  {:>4} {: >21.13}  {: >21.13}".format(i,E,dE))

			self.E = E
			self.C = C
			self.e = e
			self.D = D 

			if dE < self.conv:
				break

		return self.E


"""
spin-blocking functions: transform from spatial orbital {x_mu} basis to spin orbital basis {x_mu alpha, x_mu beta}
"""
# 1-electron integrals
def block_oei(A):
    A = np.matrix(A)
    O = np.zeros(A.shape)
    return np.bmat( [[A,O],[O,A]] )

# 2-electron integrals
def block_tei(T):
    t = np.array(T)
    n = t.shape[0]
    I2 = np.identity(2)
    T = np.zeros( (2*n,2*n,2*n,2*n) )
    for p in range(n):
        for q in range(n):
            T[p,q] = np.kron( I2, t[p,q] )
            T[n:,n:] = T[:n,:n]
    return T

if __name__ == '__main__':
	
	config = cfp.ConfigParser()
	config.read('Options.ini')

	molecule = psi4.geometry( config['DEFAULT']['molecule'] )
	molecule.update_geometry()

	basis = psi4.core.BasisSet.build(molecule, "BASIS", config['DEFAULT']['basis'],puream=0)
	mints = psi4.core.MintsHelper(basis)

	maxiter   = int( config['SCF']['maxIter'] )
	conv      = float( config['SCF']['conv']  )

	uhf = UHF(molecule,mints,maxiter,conv)
	uhf.computeEnergy()
