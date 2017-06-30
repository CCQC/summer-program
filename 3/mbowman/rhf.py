import psi4.core
import numpy as np
import scipy.linalg as spla
import math

class RHF(object):

	def __init__(self, mol, mints, convCrit = 10, maxIter = 200):
		"""
		:param mol: 
		:param mints:
		:param convCrit: criteria for converge, x corresponds to maximum difference of 10^-x
		"""
			
		self.mol = mol
		self.mints = mints
		self.convCrit = convCrit
		self.maxIter = maxIter
		self.E = 0.0

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
		self.Eold = 1.0
		self.Dold = np.zeros((self.norb,self.norb))  
		self.iter = 0
		while not convCond:
			self.I2SC()	
			if (np.absolute((self.Eold - self.E)) < 10**(-self.convCrit)  and self.iter > 1):
				print "Converged"
				print "Energy: " + str(self.E)
				print "Iterations: " + str(self.iter)	
				convCond = True
			elif self.iter >= self.maxIter:
				print "Failed to converge"
				break
		
		 
		 
	def I2SC(self):
		"""
		Iteration to(2) Self Consistency, this function is called iteratively until the energy converges
		"""
		pass
		self.Eold = self.E
		self.Dold = self.D

		#Step 1: Build Fock matrix (ie "giving a fock")
		self.J = np.einsum('prqs,rs->pq', self.g, self.D) #Sum of Columb Operators
		# print self.J
		self.K = np.einsum('prsq,rs->pq', self.g, self.D) #Sum of Exchange Operators
		# print self.K
		self.F = self.H + 2*self.J - self.K	#Fock Operator
		# print self.F	
		
		#Step 2: Transform Fock to orthogonalized AO basis
		self.FOrtho = np.einsum('ik,kj->ij', (np.einsum('ik,kj->ij',self.X,self.F)) , self.X) #Orthogonalized Fock matrix (XFX)
		# print self.FOrtho	
		
		#Step 3: Diagonilze Orthonalized Fock 
		self.orbitalEnergies, self.COrtho = np.linalg.eigh(self.FOrtho) #Orbital Energies and MO coefficients
		# print self.COrtho
		
		#Step 4: Backtransform MO Coefficients
		self.C = np.einsum('ik,kj->ij',self.X,self.COrtho) #Original basis
		#print self.C
		
		#Step 5: Rebuild density matrix
		self.Cocp = self.C[:,:self.norb]
		self.D = np.einsum('pi,qi->pq',self.Cocp,np.conj(self.Cocp)) #New density matrix 
		#print self.D
		
		#Step 6: Calculate energy (For convenience steps 5 and 6 have been rearranged wrt the inst
                HF = self.H + self.F #sum of Hamiltonian and Fock matrices
                self.E = np.einsum('pq,qp',HF,self.D) #New energy, congratulations
	
		#Step 7: Increase iteration count
		self.iter +=1	
		print "Energy: " + str(self.E)
		print "Iteration: " + str(self.iter)	

