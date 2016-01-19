from psi4_helper import get_nbf,get_nocc,get_conv,get_maxiter
import numpy as np
import scipy.linalg as la

class UHF:
    """
    Class for unrestricted (spin-orbital) Hartree-Fock computations
    """

    def __init__(self,mol,mints):
        """
        Initialize uhf class
        :param mol: a psi4 molecule object
        :param mints: a molecular integrals object (from psi4 MintsHelper)
        """

        self.nbf = get_nbf(mints)
        self.norb = 2* self.nbf
        self.nocc = get_nocc(mol)
        self.conv = get_conv()
        self.maxiter = get_maxiter()

        self.S = block_oei( mints.ao_overlap() )                      ##  S = overlap matrix
        self.X = np.matrix(la.funm(self.S,lambda x : x**(-0.5)))      ##  S^{-1/2} diagonalizes the Fock matrix

        self.T = block_oei( mints.ao_kinetic() )
        self.V = block_oei( mints.ao_potential() )
        self.Vnu = mol.nuclear_repulsion_energy()                     ##  Nuclear repulsion energy
        self.H = self.T + self.V
        
        G = block_tei(np.array( mints.ao_eri() ) )                    ##  4D ERI array (chemists' notation)
        self.g = G.transpose((0,2,1,3))-G.transpose(0,2,3,1)                                     ##  ( m n | r s )  --> < m r || n s >

        self.D = np.matrix(np.zeros(self.S.shape))                    ##  Density matrix

        self.E = 0.0

    def compute_energy(self):
        """
        print E and dE from each iteration
        return uhf energy
        """

        g, H, X = self.g, self.H, self.X

        for i in range(self.maxiter):
            v = np.einsum("mnrs,ns->mr", g, self.D)
            F = H+ v                                                  ##  build fock matrix (1st guess is F = Hcore) 
            e,tC = la.eigh(X * F * X)                                 ##  eigenvalues, vectors of transformed Fock matrix
            C = X * tC                                                ##  backtransform C
            oC = C[:,:self.nocc]                                      ##  occupied MO coefficients
            self.D = oC * oC.T                                             ##  density matrix (of MO coefficients)

            E0 = self.E
            E = np.trace( (H+0.5*v)*self.D) + self.Vnu           ##  HF energy
            self.TF = X* ( (H + 0.5*v)*self.D ) *X

            dE = np.fabs(E-E0)
            print("UHF  {:>4} {: >21.13}  {: >21.13}".format(i,E,dE))

            ##  object variables we changed in this iteration
            self.E, self.e, self.C, self.D = E, e, C, self.D

            if dE < self.conv: break

        self.Fmo = C.T * F * C
        return self.E
         

"""
spin-blocking functions: transform from spatial orbital {x_mu} basis to spin orbital basis {x_mu alpha, x_mu beta}
"""

# 1-electron integrals
def block_oei(A):
    A = np.matrix(A)
    O = np.zeros(A.shape)
    return np.bmat( [[A,O],[O,A]] )     # bmat makes block matrices !!! 


# 2-electron integrals
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
