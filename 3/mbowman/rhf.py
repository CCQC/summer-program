import psi4.core
import numpy as np
import scipy.linalg as spla

class RHF(object):

	def __init__(self, mol, mints, convCrit = 10):
		"""
		:param mol: 
		:param mints:
		:param convCrit: criteria for converge, x corresponds to maximum difference of 10^-x
		"""
			
		self.mol = mol
		self.mints = mints
		

		#Step 1: Read nuclear repulsion energy from molecule and atomic integrals from MintsHelper	
		self.VNuc = mol.nuclear_repulsion_energy() #nuclear repulsion energy
		self.S = np.array(mints.ao_overlap()) #overlap integrals
		self.T = np.array(mints.ao_kinetic()) #kinetic energy integrals
		self.V = np.array(mints.ao_potential()) #electron-nuclear attraction integrals
		self.g = np.array(mints.ao_eri()) #electron-electron repulsion integrals
		self.norb = mints.basisset().nbf() #number of orbits (defines size of arrays) 	
		#Step 2: Form orthogonalizer (X = S^-1/2)
		self.X = spla.inv(spla.sqrtm(self.S))
		
		self.H = self.T + self.V #Hamiltonian
		self.I = np.identity(self.norb)	
		#Step 3: Set D = 0 as "core" guess
		self.D = np.zeros((self.norb,self.norb))
		
		#Iteration to SC
		convCond = False #convergence condition initially set to false once the energy and density matrix converge it
		self.Eold = 0
		self.Energ = 0 
		while not convCond:
			self.I2SC()	
			if ((self.Eold - self.Energ) < 10**(-convCrit)):
			 	convCond = True
		 
		 
	def I2SC(self):
		"""
		Iteration to(2) Self Consistency, this function is called iteratively until the energy converges
		"""
		pass

		#Step 1: Build Fock matrix (ie "giving a fock")
		J = np.einsum('prqs,rs->pq', self.g, self.D) #Sum of Columb Operators
		K = np.einsum('prsq,rs->pq', self.g, self.D) #Sum of Exchange Operators
		F = H + 2*J - K	#Fock Operator
		
		#Step 2: Transform Fock to orthogonalized AO basis
		FOrtho = np.einsum('ik,kj->ij',(np.einsum('ik,kj->ij',X,F), X) #Orthogonalized Fock matrix (XFX)
			
		
