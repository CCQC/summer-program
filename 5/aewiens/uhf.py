from psi4_helper import get_nbf,get_nocc,get_conv,get_maxiter
import numpy as np
from scipy import linalg as la

class UHF:

   def __init__(self,mol,mints):
      self.nbf = get_nbf(mints)
      self.norb = 2* self.nbf
      self.nocc = get_nocc(mol)
      self.conv = get_conv()
      self.maxiter = get_maxiter()
      m = self.nbf
      N = self.norb

      # overlap matrix
      S = mints.ao_overlap()
      self.Z = block_oei(S)
      self.X = np.matrix(la.funm(self.Z,lambda x : x**(-0.5)))

      # KE matrix
      T = mints.ao_kinetic()
      self.T = block_oei(T)

      # PE matrix
      V = mints.ao_potential()
      self.V = block_oei(V)

      self.Vnu = mol.nuclear_repulsion_energy()

      # Hcore
      self.H = self.T + self.V

      # ERI matrix
      G = np.array( mints.ao_eri() )                                ####  mxmxmxm tensor
      self.G = block_tei(G)                                         ####  chemists' notation
      self.g = self.G.swapaxes(1,2) - self.G.swapaxes(1,3)          ####  ( m n | r s ) - ( m s | n r ) ->   < m r | n s >  -  < m r | s n >

      # Density matrix
      self.D = np.matrix(np.zeros(self.Z.shape))

      self.E = 0.0

   def get_energy(self):
      """
      print E and dE from each iteration
      return uhf energy
      """

      ####  rename object variables
      g, D, H, X, Vnu, E = self.g, self.D, self.H, self.X, self.Vnu, self.E

      for i in range(self.maxiter):
         v = np.einsum("mnrs,ns->mr", g, D)
         F = H+ v                                     # build fock matrix
         e,tC = la.eigh(X * F * X)                    #  eigenvalues, vectors of diagonalized Fock matrix
         C = X * tC                                   #  backtransform C
         oC = C[:,:self.nocc]                         #  occupied MO coefficients
         D = oC * oC.T                                #  density matrix (of MO coefficients)

         E0 = E
         E = np.trace( (self.H+0.5*v)*self.D) + self.Vnu        #  HF energy

         dE = np.fabs(E-E0)
         print("UHF  {:>4} {: >21.13}  {: >21.13}".format(i,E,dE))

         #### save the object variables we changed in this iteration
         self.D = D
         self.E = E
         self.e = e
         self.C = C

         if dE < self.conv: break
      print self.e
      return self.E
         

"""
spin-blocking functions: transform from spatial orbital {x_mu} basis to spin orbital basis {x_mu alpha, x_mu beta}
"""

# block one-electron integrals
def block_oei(A):
   A = np.matrix(A)
   O = np.zeros(A.shape)
   return np.bmat( [[A,O],[O,A]] )     # bmat makes block matrices !!! 

# block two-electron integrals
def block_tei(T):
   t = np.array(T)
   n = t.shape[0]
   I2 = np.identity(2)
   T = np.zeros( (2*n,2*n,2*n,2*n) )
   for p in range(n):
      for q in range(n):
         T[p,q] = np.kron( I2, t[p,q] )
         T[n:,n:] = T[:n,:n]
   return T
