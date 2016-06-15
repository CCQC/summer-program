import numpy as np
from psi4_helper import get_docc, get_nbf, get_conv, get_maxiter

class RHF:
# RHF class for obtaining the HF energy of a closed shell system

    def __init__(self,mol,mints):
        """
        initialize rhf class
        takes objects:
        1. mol: a psi4 molecule, get_active_molecule()
        2. mints: molecular integrals from libmints in psi4
        """
        self.docc = get_docc(mol)
        self.nbf = get_nbf(mints)
        self.conv = get_conv()
        self.maxiter = get_maxiter()

        S = np.matrix(mints.ao_overlap() )
        self.X = np.matrix(np.linalg.funm(S,lambda x : x**(-0.5)))

        T = np.matrix(mints.ao_kinetic() )
        V = np.matrix(mints.ao_potential() )
        self.H = T+V
        self.Vnu = mol.nuclear_repulsion_energy()

        self.G = np.array(mints.ao_eri() ).transpose((0,2,1,3))
        self.D = np.zeros((self.nbf,self.nbf))
        self.E = 0

    def get_energy(self):
        """
        print E and dE from each iteration
        return rhf energy
        """

        docc, nbf, conv, maxiter = self.docc, self.nbf, self.conv, self.maxiter
        X, H, G, Vnu, D, E = self.X, self.H, self.G, self.Vnu, self.D, self.E

        for i in range(maxiter):
            J= np.einsum("ikjl,kl->ij",G,D)
            K = np.einsum("iklj,kl->ij",G,D)
            F = H+J-0.5*K                           ####  build fock matrix
            tF = X*F*X                              ####  diagonalize fock matrix
            e, tC = np.linalg.eigh(tF)              ####  eigenvalues & eigenvectors of fock matrix
            C = X*tC                                ####  back-transform the eigenvectors
            oC = C[:,:docc]                         ####  save occupied block of C matrix
            D = 2*oC*oC.T                           ####  build density matrix

            E0 = E
            E = np.trace(0.5*(H+F)*D)+Vnu
            dE = np.fabs(E-E0)

            print( "{0:20.12f}{1:20.12f}".format(E,dE) )
            if dE < self.conv : break

            #### save object variables that changed during this iteration
            self.E, self.e, self.C, self.D = E, e, C, D

        return self.E
