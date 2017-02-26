#!/usr/bin/env python
import configparser, psi4, numpy as np

class UHF:

	def __init__(self, options):

		mol   = psi4.geometry( options['DEFAULT']['molecule'] )
		mol.update_geometry()

		self.basisName = options['DEFAULT']['basis']
		self.basis     = psi4.core.BasisSet.build(mol, "BASIS", self.basisName, puream=0)
		self.getIntegrals() 

		mult  = mol.multiplicity()
		nelec = self.getNelec(mol)
		self.Vnu = mol.nuclear_repulsion_energy()

		self.maxiter = int( options['SCF']['max_iter'] )
		self.conv    = 10**( -int(options['SCF']['conv']) )
		self.dConv   = 10**( -int(options['SCF']['d_conv']) )

		self.Na  = int( 0.5*(nelec+mult-1) )
		self.Nb  = nelec - self.Na
		self.Da  = np.zeros((self.X.shape))
		self.Db  = np.zeros((self.X.shape))
		self.E   = 0.0
		
		self.converged = False
		self.options   = options


	def getNelec(self,mol):

		char = mol.molecular_charge()
		nelec = -char
		for A in range(mol.natom()):
			nelec += mol.Z(A)

		return int(nelec)


	def getIntegrals(self):

		mints    = psi4.core.MintsHelper( self.basis )
		self.V   = np.array( mints.ao_potential() )
		self.T   = np.array( mints.ao_kinetic() )
		self.S   = np.array( mints.ao_overlap() )

		G = np.array( mints.ao_eri() )
		self.G = G.transpose((0,2,1,3))

		S = mints.ao_overlap()
		S.power( -0.5, 1.e-16 )
		self.X = np.array( S )


	def computeEnergy(self):

		H = self.T + self.V
		G = self.G
		Da = self.Da
		Db = self.Db
		X  = self.X

		print('\n      Iter         Energy                 ΔE                   ‖ΔD‖')
		print('----------------------------------------------------------------------------')

		for i in range( self.maxiter ):

			va = np.einsum("mnrs,ns->mr",G,Da+Db) - np.einsum("mnsr,ns->mr",G,Da)
			vb = np.einsum("mnrs,ns->mr",G,Da+Db) - np.einsum("mnsr,ns->mr",G,Db)
			Fa = H + va
			Fb = H + vb

			tFa = X@Fa@X
			tFb = X@Fb@X

			ea, tCa = np.linalg.eigh(tFa)
			eb, tCb = np.linalg.eigh(tFb)

			Ca = X@tCa
			Cb = X@tCb

			oCa = Ca[:,:self.Na]
			oCb = Cb[:,:self.Nb]
			Da  = oCa@oCa.T
			Db  = oCb@oCb.T

			E   = np.trace( (H+0.5*va)@Da ) + np.trace( (H+0.5*vb)@Db ) + self.Vnu
			dE  = np.fabs( E - self.E )
			dDa = np.fabs( np.linalg.norm(Da) - np.linalg.norm(self.Da) )
			dDb = np.fabs( np.linalg.norm(Db) - np.linalg.norm(self.Db) )
			dD  = (dDa+dDb)/2

			if dE < self.conv and dD < self.dConv:
				self.converged = True
				break

			print("UHF {:>4d}{: >21.12f}{: >21.12f}{: >21.12f}".format( i, E, dE, dD ))

			self.Da = Da
			self.Db = Db
			self.E  = E

		self.writeOutput() 
		return self.E


	def psiSCF(self):

		psi4.core.set_output_file("output.dat",False)

		psi4.set_options({'basis': self.options['DEFAULT']['basis'],
						  'scf_type': 'pk',
						  'reference': 'uhf',
						  'puream': 0,
						  'print': 0 })

		return psi4.energy('scf')


	def writeOutput(self):

		print('\n----------------------------------------------------------------------------')
		if self.converged:
			print("UHF procedure has converged.")
			print("Final UHF energy:{:20.11f}".format(self.E) )
			print("Compare to Psi4: {:20.11f}".format( self.psiSCF() ) )

		print('----------------------------------------------------------------------------')
		if not self.converged:
			print("UHF procedure did not converge. Sorry")
		


if __name__ == '__main__':

	config = configparser.ConfigParser()
	config.read('Options.ini')
	uhf   = UHF(config)
	uhf.computeEnergy()
