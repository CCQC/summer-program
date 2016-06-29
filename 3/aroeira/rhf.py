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

        self.g = np.array(mints.ao_eri())

        # Determine the number of electrons and the number of doubly occupied orbitals
        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        if mol.multiplicity() != 1 or self.nelec % 2:
            raise Exception("This code only allows closed-shell molecules")
        self.ndocc = self.nelec / 2

        self.maxiter = psi4.get_global_option('MAXITER')
        self.e_convergence = psi4.get_global_option('E_CONVERGENCE')
        self.nbf = mints.basisset().nbf()
        self.vu = np.matrix(np.zeros((self.nbf, self.nbf)))
    
    def build_vu(self, D):
        """
        Construct a v matrix using the density matrix and the four index integrals
        """

        r = range(int(self.nbf))
        self.vu = np.matrix(np.zeros((self.nbf, self.nbf)))
        for u in r:
            for v in r:
                for p in r:
                    for q in r:
                        self.vu[u, v] += (2 * self.g[u, v, p, q] - self.g[u, p, v, q]) * D[q, p] 

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
            self.build_vu(D)
            E1 = np.sum((2 * np.array(h) + np.array(self.vu))*np.array(D.T)) + self.V_nuc
            psi4.print_out('Iteration {:<d}   {:.10f}    {:.10f}\n'.format(count, E1, E1-E0))
            if abs(E1 - E0) < self.e_convergence:
                psi4.print_out('\nFinal Energy: {:<5.10f}'.format(E1))
                break
            else:
                E0 = E1
        else:
            psi4.print_out('\n:(   Does not converge   :(')
