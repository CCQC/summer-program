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
        nbf = get_nbf(mints)
        self.norb = 2* nbf
        self.nocc = get_nocc(mol)
        self.conv = get_conv()
        self.maxiter = get_maxiter()

        T = block_oei( mints.ao_kinetic() )
        V = block_oei( mints.ao_potential() )
        S = block_oei( mints.ao_overlap() )                      
        G = block_tei(np.array( mints.ao_eri() ) )                    ##  4D ERI array (chemists' notation)

        self.H = T + V
        self.X = np.matrix(la.funm(S,lambda x : x**(-0.5)))           ##  S^{-1/2} diagonalizes the Fock matrix
        self.Vnu = mol.nuclear_repulsion_energy()                     ##  Nuclear repulsion energy
        self.g = G.transpose((0,2,1,3))-G.transpose((0,2,3,1))        ##  ( m n | r s )  --> < m r || n s >
        self.D = np.matrix(np.zeros( S.shape ))                    

        self.E = 0.0

    def compute_energy(self):
        """
        print E and dE from each iteration
        return uhf energy
        """

        g, H, X, D = self.g, self.H, self.X, self.D

        for i in range(self.maxiter):
            v = np.einsum("mnrs,ns->mr", g, self.D)
            F = H+ v                                                  ##  build fock matrix (1st guess is F = Hcore) 
            e,tC = la.eigh(X * F * X)                                 ##  eigenvalues, vectors of transformed Fock matrix
            C = X * tC                                                ##  backtransform C
            oC = C[:,:self.nocc]                                      ##  occupied MO coefficients
            D = oC * oC.T                                             ##  density matrix (of MO coefficients)

            E0 = self.E
            E = np.trace( (H+0.5*v)*D) + self.Vnu           ##  HF energy
            dE = np.fabs(E-E0)
            print("UHF  {:>4} {: >21.13}  {: >21.13}".format(i,E,dE))

            ##  object variables we changed in this iteration
            self.E, self.e, self.C, self.D = E, e, C, D

            if dE < self.conv: break

        return self.E
         

"""
spin-blocking functions: transform from spatial orbital {x_mu} basis to spin orbital basis {x_mu alpha, x_mu beta}
"""
# 1-electron integrals
def block_oei(A):
    A = np.matrix(A)
    O = np.zeros(A.shape)
    return np.bmat( [[A,O],[O,A]] )

# 2-electron integrals
# what on earth is happening here
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
