import psi4.core
import numpy as np
import scipy.linalg as spla

class RHF(object):

	def __init__(self, mol, mints):
		"""
		:param mol: 
		:param mints:
		"""
			
		self.mol = mol
		self.mints = mints
		

		#Step 1 Read nuclear repulsion energy from molecule and atomic integrals from MintsHelper	
		self.VNuc = mol.nuclear_repulsion_energy() #nuclear repulsion energy
		self.S = np.array(mints.ao_overlap()) #overlap integrals
		self.T = np.array(mints.ao_kinetic()) #kinetic energy integrals
		self.V = np.array(mints.ao_potential()) #electron-nuclear attraction integrals
		self.g = np.array(mints.ao_eri()) #electron-electron repulsion integrals

		#Step 2 Form orthogonalizer (X = S^-1/2)
		self.X = spla.inv(spla.sqrtm(self.S))
		


