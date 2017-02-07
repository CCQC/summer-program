#!/usr/bin/env python3

import psi4, configparser, numpy as np

class RHF(object):

	def __init__(self,molecule,mints):

		self.docc = self.getDocc(molecule)       
		self.Vnu  = molecule.nuclear_repulsion_energy()
		self.E    = 0

		self.getIntegrals(mints)
		self.D = np.zeros_like( self.S )

		self.converged = False


	def getDocc(self,mol):
		char = mol.molecular_charge()
		nelec = -char
		for A in range(mol.natom()):
			nelec += mol.Z(A)
		return int(nelec//2)


	def getIntegrals(self,mints):

		self.V = np.array( mints.ao_potential() )
		self.T = np.array( mints.ao_kinetic() )
		self.G = np.array( mints.ao_eri() ).transpose((0,2,1,3))

		S = mints.ao_overlap()
		S.power(-0.5, 1.e-16)
		self.X = S.to_array() 
		self.S = np.array( mints.ao_overlap() )


	def computeEnergy(self):

		H  = self.T + self.V
		G  = self.G
		X  = self.X
		S  = self.S
		D  = self.D
		docc = self.docc

		for i in range(100):

			J= np.einsum("ikjl,kl",G,D)
			K = np.einsum("iklj,kl",G,D)
			F = H+2*J-K

			tF = X.dot(F.dot(X))

			e, tC = np.linalg.eigh(tF)

			C    = X.dot(tC)
			Cocc = C[:,:docc]
			D    = Cocc.dot(Cocc.T)

		
			E0 = self.E
			E  = np.trace((H+F) @ D)+ self.Vnu
			dE = np.fabs(E-E0)

			if __name__ == '__main__':
				print("RHF  {:>4} {: >21.11}  {: >21.11}".format(i,E,dE))

			if dE < 1e-12:
				self.converged = True
				self.writeOutput()
				break

			self.D = D
			self.E = E

		return E


	def writeOutput(self):

		if __name__ == '__main__':

			if self.converged:
				print("\nRHF procedure has converged.\nFinal RHF energy:{:20.11f}".format(self.E) )
			if not self.converged:
				print("RHF procedure did not converge. Sorry")




if __name__ == '__main__':

	config = configparser.ConfigParser()
	config.read('Options.ini')
	molecule   = psi4.geometry( config['DEFAULT']['molecule'] )
	molecule.update_geometry()
	basis = psi4.core.BasisSet.build(molecule, "BASIS", config['DEFAULT']['basis'],puream=0)
	mints = psi4.core.MintsHelper(basis)
	rhf   = RHF(molecule,mints)
	rhf.computeEnergy()
