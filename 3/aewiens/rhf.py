import numpy as np
from scipy import linalg as la
from psi4_helper import get_docc, get_nbf, get_conv, get_maxiter

class RHF:
# RHF class for obtaining the HF energy of a closed shell system

    def __init__(self,mol,mints):
        """
        initialize the class
        takes objects:
        1. mol: a psi4 molecule, get_active_molecule()
        2. mints: molecular integrals from libmints in psi4
        """
        self.docc = get_docc(mol)
        self.nbf = get_nbf(mints)
        self.conv = get_conv()
        self.maxiter = get_maxiter()
        self.S = np.matrix(mints.ao_overlap() )
        self.X = np.matrix(la.funm(self.S,lambda x : x**(-0.5)))
        self.T = np.matrix(mints.ao_kinetic() )
        self.V = np.matrix(mints.ao_potential() )
        self.H = self.T+self.V
        self.G = np.array(mints.ao_eri() ).swapaxes(1,2)
        self.Vnu = mol.nuclear_repulsion_energy()
        self.D = np.matrix(np.zeros(self.S.shape))
        self.E = 0

    def get_energy(self):
        """
        print E and dE from each iteration
        return rhf energy
        """

        #rename class variables
        docc, nbf, conv, maxiter = self.docc, self.nbf, self.conv, self.maxiter
        S, X, T, V, H, G, Vnu, D = self.S, self.X, self.T, self.V, self.H, self.G, self.Vnu, self.D

        E = self.E

        for i in range(maxiter):
            J= np.einsum("ijkl,jl->ik",G,D)
            K = np.einsum("ijkl,jk->il",G,D)
            F = H+J-0.5*K      # start with F = Hcore
            tF = X*F*X
            e, tC = la.eigh(tF)
            C = X*tC
            oC = C[:,:docc]
            D = 2*oC*oC.T
            E0 = E
            E = np.trace(0.5*(H+F)*D)+Vnu
            dE = np.fabs(E-E0)
            self.E = E
            print("%20.12f%20.12f" % (E, dE))
            if dE < self.conv : break
        return self.E
