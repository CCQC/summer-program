import psi4
import numpy as np
import scipy.linalg as la


class RHF:
    """
    Restricted Hartree-Fock class for obtaining the restricted Hartree-Fock
    energy
    """

    def __init__(self, mol, mints):
        """
        Initialize the rhf
        :param mol: a Psi4 molecule object
        :param mints: a molecular integrals object (from MintsHelper)
        """
        self.mol = mol
        self.mints = mints

        self.V_nuc = mol.nuclear_repulsion_energy()
        self.T = np.matrix(mints.ao_kinetic())
        self.S = np.matrix(mints.ao_overlap())
        self.V = np.matrix(mints.ao_potential())

        self.g = np.array(mints.ao_eri()).swapaxes(1,2)

        # Determine the number of electrons and the number of doubly occupied orbitals
        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        if mol.multiplicity() != 1 or self.nelec % 2:
            raise Exception("This code only allows closed-shell molecules")
        self.ndocc = int(self.nelec / 2)

        self.maxiter = psi4.core.get_global_option('MAXITER')
        self.e_convergence = psi4.core.get_global_option('E_CONVERGENCE')
        self.nbf = mints.basisset().nbf()
        self.vu = np.matrix(np.zeros((self.nbf, self.nbf)))
    
    def compute_energy(self):
        """
        Compute the rhf energy
        :return: energy
        """
        X = np.matrix(la.inv(la.sqrtm(self.S)))
        D = np.matrix(np.zeros((self.nbf, self.nbf)))
        h = self.T + self.V
        E0 = 0
        for count in range(self.maxiter):
            F = h + self.vu
            Ft = X * F * X
            e, Ct = la.eigh(Ft)
            C = X * np.matrix(Ct)
            DOC = np.matrix(C[:,:self.ndocc])
            D = DOC*DOC.T
            G = 2*self.g - self.g.swapaxes(2,3)
            self.vu = np.einsum('upvq,qp->uv', G, D) 
            E1 = np.sum((2 * np.array(h) + np.array(self.vu))*np.array(D.T)) + self.V_nuc
            psi4.core.print_out('Iteration {:<d}   {:.10f}    {:.10f}\n'.format(count, E1, E1-E0))
            if abs(E1 - E0) < self.e_convergence:
                psi4.core.print_out('\nFinal HF Energy: {:<5.10f}'.format(E1))
                self.C = C
                self.epsi = e
                self.ehf = E1
                break
            else:
                E0 = E1
        else:
            psi4.core.print_out('\n:(   Does not converge   :(')

class MP2(RHF):

    def mp2_energy(self):
        E2 = 0
        self.G = np.einsum('uvpq,uP,vQ,pR,qS->PQRS', self.g, self.C, self.C, self.C, self.C)
        for I in range(self.ndocc):
            for J in range(self.ndocc):
                for A in range(self.ndocc, len(self.epsi)):
                    for B in range(self.ndocc, len(self.epsi)):
                        E2 += (self.G[I, J, A, B]*(2*self.G[I, J, A, B] - self.G[I, J, B, A]))/(self.epsi[I] + self.epsi[J] - self.epsi[A] - self.epsi[B])
        self.emp2 = E2
        psi4.core.print_out('\nMP2 Energy: {:<5.10f}'.format(self.emp2))
        psi4.core.print_out('\nTotal Energy: {:<5.10f}'.format(self.ehf + self.emp2))
        return self.emp2
