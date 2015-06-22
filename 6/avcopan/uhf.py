import psi4
import numpy as np
import scipy.linalg as la

class UHF:

    def __init__(self, mol, mints):
        self.mol = mol
        self.mints = mints

        self.V_nuc = mol.nuclear_repulsion_energy()
        self.T = block_matrix_aa(mints.ao_kinetic())
        self.S = block_matrix_aa(mints.ao_overlap())
        self.V = block_matrix_aa(mints.ao_potential())

        self.g = block_4darray(mints.ao_eri())

        print self.S.round(2)


def block_matrix_aa(A):
  A  = np.matrix(A)
  I2 = np.identity(2)
  return np.matrix( np.kron(I2, A) )


def block_4darray(T):
  t  = np.array(T)
  n  = t.shape[0]
  I2 = np.identity(2)
  T  = np.zeros((2*n,2*n,2*n,2*n))
  for p in range(n):
    for q in range(n):
      T[p,q] = np.kron(I2, t[p,q])
  T[n:,n:] = T[:n,:n]
  return T

