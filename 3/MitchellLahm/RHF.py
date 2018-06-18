import numpy as np
import scipy.linalg as la
import psi4.core

class RHF():
	def __init__(self, mol, mints):
		# Reading in the molecule, molec ints, and Nuclear repulsion energy
		self.mol = mol
		self.mints = mints
		self.E = 0
		self.Vnu = mol.nuclear_repulsion_energy()
		# Reading in other pertinent data
		self.natom = mol.natom()
		self.charge = mol.molecular_charge()
		self.norb = mints.basisset().nbf()
		# Defining molecular integral matrices
		self.S = np.array(mints.ao_overlap())    # overlap integrals
		self.T = np.array(mints.ao_kinetic())    # kinetic energy integrals
		self.V = np.array(mints.ao_potential())  # electron-nuclear attraction integrals
		self.g = np.array(mints.ao_eri())        # electron-electron repulsion integrals
		self.g = self.g.transpose((0,2,1,3))     # Converting to physicist's notation
		self.X = la.inv(la.sqrtm(self.S))        # Generate orthogonalizer
		self.D = np.zeros((self.norb,self.norb)) # Initialize empty density matrix
	# Iteration to Self-Consistency
	def compute_energy(self):
		# Defining necessary terms
		dE = 1
		count = 0
		nelec = -self.charge + sum(self.mol.Z(A) for A in range(self.natom))
		# any open shell systems will be forced to occupy an integer number of orbitals
		nocc = int(nelec / 2.0)
		X = self.X
		D = self.D
		Vnu = self.Vnu
		h = self.T + self.V
		j = self.g
		k = self.g.transpose((0,1,3,2))
		while(dE > 1e-11 and count < 50):
			# Building the fock matrix
			u1 = np.tensordot(j, D, axes = [(1,3),(1,0)])
			u2 = -0.5 * np.tensordot(k, D, axes = [(1,3),(1,0)])
			u = u1 + u2
			f = h + u
			# Transforming f to orthogonalized AO basis
			f = X.dot(f).dot(X)
			# Grabbing orbital energies, MO coefficients
			e,C = la.eigh(f)
			# Transforming ~C back to C
			C = X.dot(C)
			self.C = C
			# Building Density matrix
			Cs = C[:,:nocc]
			D = 2 * Cs.dot(Cs.T)
			# Evaluate energy
			Ee = 0.5 * u + h
			#Ee = (Ee * D.T).sum()
			Ee = np.trace(Ee.dot(D))
			E = Ee + Vnu
			# Checking convergence condition
			dE = abs(E-self.E)
			self.E = E
			count += 1
			psi4.core.print_out('Iteration: %2.0f dE: %17.13f RHF-Energy: %16.10f \n'%(count,dE,E))
		self.e,self.C = np.linalg.eigh(f)
		

		
		
