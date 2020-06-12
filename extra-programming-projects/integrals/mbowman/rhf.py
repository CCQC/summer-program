import numpy as np
import scipy.linalg as spla
import math
import masses
import molecule
import integrals
	
np.set_printoptions(threshold=np.inf, precision=3, linewidth=200, suppress=True)

class RHF(object):

	def __init__(self, mol_str,convCrit=10,maxIter=100):
		"""
		Initializes molecule, molecular properties, molecular integrals,
		before entering a loop to calculate the RHF energy 

		:param mol: molecule object, specifies geometry, charge, and multiplicity 
		:param mints: molecular integral helper, generates various molecular intergrals
		:param convCrit: criteria for converge, x in input corresponds to maximum difference of 10^-x
		:param maxIter: maximum number of iterations to obtain self consistence
		"""
			
		mol = molecule.Molecule(mol_str) 
		#Step 1: Read nuclear repulsion energy from molecule and atomic integrals from MintsHelper	
		mints = integrals.Integral(mol_str)	

		self.VNuc = mints.VNuc 		#nuclear repulsion energy
		self.S = np.array(mints.S) 			#overlap integrals
		self.T = np.array(mints.T) 			#kinetic energy integrals
		self.V = np.array(mints.V) 		#electron-nuclear attraction integrals
		self.g = np.array(mints.G).transpose(0,2,1,3) 	#electron-electron repulsion integrals, transposed from (pq|rs) to <pr|qs>
		#The transpose is very important, because from here on forward physicist notation is assumed!
		
		#Step 1.5: Calculate Hamiltonian and orbital information 
		self.H = self.T + self.V 		#Hamiltonian
		self.E = 0.0 				#RHF energy
		self.norb = mints.K 	#number of orbits (defines size of arrays) 	
		nelec = mol.molecular_charge + sum(mol.charges) 	#number of electrons
		self.nocc = int(nelec / 2)		#number of occupied orbitals
					
		#Step 2: Form orthogonalizer (X = S^-1/2)
		self.X = np.matrix(spla.inv(spla.sqrtm(self.S)))	#Orthogonalizer S^-1/2
		
		#Step 3: Set D = 0 as "core" guess
		self.D = np.zeros((self.norb,self.norb))		#Density Matrix
		
		#Iteration to SC
		convCond = False 				#convergence condition 
		self.Eold = 1.0					#Previous calculated energy
		self.Dold = np.zeros((self.norb,self.norb))  	#Previous density matrix
		self.iter = 0					#Iteration count
		line = "+----+---------------------+------------+"
		print line
		print  "|iter| Energy              | dE         |"
		print line
		while not convCond:
			self.I2SC() #Steps 1-7
			#Check convergence, if the energy and density matrix are both within a threshold of one another, then it has converged
			#Additionally the program must iterate twice to avoid the condition when all four variables are initially null
			if (np.absolute((self.Eold - self.E)) < 10.0**(-convCrit) and np.absolute(spla.norm(self.D) - spla.norm(self.Dold)) < 10.0**(-convCrit) and self.iter > 1):
				print line
				print  "| Converged                             |"
				convCond = True
			elif self.iter >= maxIter:
				print line
				print  "| Failed to converge                    |"
				break
		print line
		 
		 
	def I2SC(self):
		"""
		Iteration to(2) Self Consistency, this function is called iteratively until the energy converges
		"""
		#Step 0: set the old values for E and D to the result of the previous iteration
		self.Eold = self.E
		self.Dold = self.D

		#Step 1: Build Fock matrix 
		self.J = np.einsum('prqs,sr->pq', self.g, self.D)	#Sum of Columb Operators
		self.K = np.einsum('prsq,sr->pq', self.g, self.D) 	#Sum of Exchange Operators
		self.F = self.H + self.J - 0.5 *self.K			#Fock Operator
		
		#Step 2: Transform Fock to orthogonalized AO basis
		self.FOrtho = np.einsum('ik,kj->ij', (np.einsum('ik,kj->ij',self.X,self.F)) , self.X) #Orthogonalized Fock matrix (XFX)
		
		#Step 3: Diagonilze Orthonalized Fock 
		self.orbitalEnergies, self.COrtho = np.linalg.eigh(self.FOrtho) #Orbital Energies and MO coefficients
		
		#Step 4: Backtransform MO Coefficients
		self.C = np.einsum('ik,kj->ij',self.X,self.COrtho) 	#Original basis
		
		#Step 5: Rebuild density matrix
		self.Cocp = self.C[:,:self.nocc]				#Truncated basis to only include occupied orbitals
		self.D = 2*np.einsum('pi,qi->pq',self.Cocp,np.conj(self.Cocp)) 	#New density matrix 
		
		#Step 6: Calculate energy 
                sum = self.H + 0.5* self.J+ 0.25*self.K 		#sum of H + 1/2 v
                self.E = np.einsum('pq,qp',sum,self.D) + self.VNuc 	#New energy, congratulations
	
		#Step 7: Print and increase iteration count
		print '|{0:>4d}| {1:>17.14f} | {2:>4.4e} |'.format(self.iter, self.E,np.absolute(self.E - self.Eold))
		self.iter +=1	

if __name__ == "__main__":
	 with open("methyl.xyz") as moleOpen:
	        ms = moleOpen.read()	
		f = RHF(ms)
