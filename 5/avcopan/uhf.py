import psi4
import numpy as np
import scipy.linalg as la


class UHF:

    def __init__(self, mol, mints):

      # compute and process integrals (blocking functions defined below)
      T = block_oei(mints.ao_kinetic())
      S = block_oei(mints.ao_overlap())
      V = block_oei(mints.ao_potential())
      G = block_tei(mints.ao_eri()) # two-electron integrals (chemist's notation!)

      # object attributes
      self.h = T + V                                       # core hamiltonian (h = T + V)
      self.g = G.transpose(0,2,1,3) - G.transpose(0,2,3,1) # antisymmetrized TEIs <mu nu || rh si> (physicist's notation)
      self.D = np.zeros(S.shape)                           # empty density matrix (core guess D=0)
      self.X = np.matrix(la.inv(la.sqrtm(S)))              # orthogonalizer (X=S^-1/2)

      self.Vnu  = mol.nuclear_repulsion_energy()
      self.nocc = int( sum(mol.Z(A) for A in range(mol.natom())) - mol.molecular_charge() ) # num e- = sum_A Z_A - mol. charge


    def compute_energy(self):
      # copy over object attributes to avoid having to write "self." a lot
      h, g, D, X, Vnu, nocc = self.h, self.g, self.D, self.X, self.Vnu, self.nocc

      self.E = 0.0

      for i in range(psi4.get_option('SCF', 'MAXITER')):
        v = np.einsum('mnrs,ns', g, D)   # e- field  v_mu,nu = sum_rh,si <mu nu||rh si> D_nu,si
        F = h + v                        # build fock matrix

        tF = X*F*X                       # transform to orthogonalized AO basis
        e, tC = la.eigh(tF)              # diagonalize
        C = X * tC                       # backtransform
        oC = C[:,:nocc]                  # get occupied MO coefficients
        D = oC * oC.T                    # build density matrix

        E  = np.trace((h + v/2)*D) + Vnu # compute energy by tracing with density matrix

        dE = E - self.E
        self.E, self.C = E, C # save MO coefficients and energy

        print('{:-3d}{:20.15f}{:20.15f}'.format(i, E, dE))
        if(np.fabs(dE) < psi4.get_global_option('E_CONVERGENCE')): break

      return self.E


# spin blocking functions: transform from spatial orbital {x_mu} basis to spin-orbital {x_mu alpha, x_mu beta} basis
# block one-electron integrals
def block_oei(A):
  A  = np.matrix(A)
  I2 = np.identity(2)
  return np.matrix( np.kron(I2, A) )

# block two-electron integrals [must be in chemist's notation, (mu nu|rh si)]
def block_tei(T):
  t  = np.array(T)
  n  = t.shape[0]
  I2 = np.identity(2)
  T  = np.zeros((2*n,2*n,2*n,2*n))
  for p in range(n):
    for q in range(n):
      T[p,q] = np.kron(I2, t[p,q])
  T[n:,n:] = T[:n,:n]
  return T

