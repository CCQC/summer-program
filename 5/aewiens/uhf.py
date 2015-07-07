from psi4_helper import get_nbf
import numpy as np

class UHF:

	def __init__(self,mol,mints):
		self.nbf = get_nbf(mints)
		self.norb = 2* self.nbf
		m = self.nbf
		N = self.norb

		# overlap matrix
		S = mints.ao_overlap()
		self.Z = block_oei(S)

		# KE matrix
		T = mints.ao_kinetic()
		self.T = block_oei(T)

		# PE matrix
		V = mints.ao_potential()
		self.V = block_oei(V)

		# Hcore
		self.H = self.T + self.V

		# ERI matrix
		G = np.array( mints.ao_eri() )    # mxmxmxm 4D tensor
		self.G = block_tei(G)
		self.g = self.G.swapaxes(1,2)   # ( m n | r s ) ->   < m r | n s > 


# spin-blocking functions: transform from spatial orbital {x_mu} basis to spin orbital basis {x_mu \alpha, x_mu \beta}

# block one-electron integrals
def block_oei(A):
	A = np.matrix(A)
	O = np.zeros(A.shape)
	return np.bmat( [[A,O],[O,A]] )     # bmat makes block matrices !!! 

# block two-electron integrals
def block_tei(T):
	t = np.array(T)
	n = t.shape[0]
	print n
	I2 = np.identity(2)
	print I2
	T = np.zeros( (2*n,2*n,2*n,2*n) )
	for p in range(n):
		for q in range(n):
			T[p,q] = np.kron( I2, t[p,q] )
	T[n:,n:] = T[:n,:n]
	return T
