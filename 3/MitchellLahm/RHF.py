import numpy as np
import scipy.linalg as la
import psi4.core
import sys


class RHF():
	def __init__(self, mol, mints, basis):
		print(str(basis))
		self.mol = mol
		self.mints = mints
		self.E = 0
		self.Vnu = mol.nuclear_repulsion_energy()
		self.natom = mol.natom()
		self.charge = mol.molecular_charge()
		self.norb = mints.basisset().nbf()
		self.S = np.array(mints.ao_overlap())    # overlap integrals
		self.T = np.array(mints.ao_kinetic())    # kinetic energy integrals
		self.V = np.array(mints.ao_potential())  # electron-nuclear attraction integrals
		self.g = np.array(mints.ao_eri())        # electron-electron repulsion integrals
		self.X = la.inv(la.sqrtm(self.S))        # Generate orthogonalizer
		self.D = np.zeros((self.norb,self.norb))
	def compute_energy(self):
		diff = 1
		count = 0
		while(diff > 1e-10 and count < 50):
			print(count)
			print(self.E)
			print(diff)
			# Building the fock matrix
			h = self.T + self.V
			j = self.g
			k = -0.5 * self.g
			u = np.einsum('ijkl,ijlk->ik',j,k).dot(self.D)
			f = h + u
			# Transforming f to orthogonalized AO basis
			f = self.X.dot(f).dot(self.X)
			# Grabbing orbital energies, MO coefficients
			self.e,self.C = np.linalg.eigh(f)
			# Transforming ~C back to C
			self.C = self.X.dot(self.C)
			# Building Density matrix
			Cs = np.matrix.conjugate(self.C)
			self.D = np.einsum('mi,vi->mv',self.C,Cs)
			# Evaluate energy
			sum = 0
			for i in range(self.norb):
				for j in range(self.norb):
					sum += (h[i][j] + 0.5 * u[i][j]) * self.D[j][i]
			E = sum + self.Vnu
			# Checking convergence condition
			diff = abs(E-self.E)
			self.E = E
			count += 1
		print(count)
		print(self.E)
		print(diff)
		print(self.e)


