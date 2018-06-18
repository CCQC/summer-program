#!/usr/bin/env python3
import psi4, configparser, numpy as np

class RHF(object):

	def __init__(self,options):

		self.mol   = psi4.geometry( options['DEFAULT']['molecule'] )
		self.mol.update_geometry()
		self.basisName = options['DEFAULT']['basis']

		self.basis = psi4.core.BasisSet.build(self.mol, "BASIS", self.basisName ,puream=0)
		self.docc  = self.getDocc(self.mol)       
		self.Vnu   = self.mol.nuclear_repulsion_energy()
		self.E     = 0

		self.maxiter = int( options['SCF']['max_iter'] )
		self.conv    = 10**( -int( options['SCF']['conv'] ))
		self.dConv   = 10**( -int( options['SCF']['d_conv'] ))

		self.getIntegrals(self.basis)
		self.D = np.zeros_like( self.S )

		self.converged = False
		self.options   = options


	def getDocc(self,mol):

		char = mol.molecular_charge()
		nelec = -char
		for A in range(mol.natom()):
			nelec += mol.Z(A)

		return int(nelec//2)


	def getIntegrals(self,basis):

		mints  = psi4.core.MintsHelper(basis)
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

		print('\n           Iter         Energy                 ΔE                    ‖ΔD‖')
		print('--------------------------------------------------------------------------------')

		for i in range(self.maxiter):

			J  = np.einsum("ikjl,kl",G,D)
			K  = np.einsum("iklj,kl",G,D)
			F  = H+2*J-K
			tF = X@F@X

			e,tC = np.linalg.eigh(tF)

			C    = X@tC
			Cocc = C[:,:docc]
			D    = Cocc@Cocc.T
		
			E0   = self.E
			E    = np.trace((H+F) @ D)+ self.Vnu
			dE   = np.fabs(E-E0)
			dD   = np.fabs( np.linalg.norm(D) - np.linalg.norm(self.D) )

			if __name__ == '__main__':
				print("@RHF-iter {:>4}{: >20.11f}{: >20.11f}{:20.11f}".format(i,E,dE,dD))

			if dE < self.conv and dD < self.dConv:
				self.converged = True
				self.writeOutput()
				break

			self.D = D
			self.E = E
			self.e = e
			self.C = C

		return E


	def psiSCF(self):

		psi4.core.set_output_file("output.dat",False)

		psi4.set_options({'basis': self.options['DEFAULT']['basis'],
						  'scf_type': 'pk',
						  'reference': 'rhf',
						  'puream': 0,
						  'print': 0 })

		return psi4.energy('scf')


	def writeOutput(self):

		if __name__ == '__main__':
			
			print('\n--------------------------------------------------------------------------------')
			if self.converged:
				print("RHF procedure has converged.")
				print("Final RHF energy:{:20.11f}".format(self.E) )
				print("Compare to Psi4: {:20.11f}".format( self.psiSCF() ) )
			print('--------------------------------------------------------------------------------')
			if not self.converged:
				print("RHF procedure did not converge. Sorry")


if __name__ == '__main__':

	config = configparser.ConfigParser()
	config.read('Options.ini')
	rhf   = RHF(config)
	rhf.computeEnergy()
