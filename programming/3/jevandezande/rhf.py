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
        :param mol: a psi4 molecule object
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

    def compute_energy(self):
        """
        Compute the rhf energy
        :return: energy
        """
        V_nuc, T, S, V, g = self.V_nuc, self.T, self.S, self.V, self.g
        nbf, ndocc = self.nbf, self.ndocc

        D = np.zeros((self.nbf, self.nbf))
        C = np.zeros((self.nbf, self.nbf))

        # S^{-1/2}
        X = la.inv(la.sqrtm(S)).view(np.matrix)

        psi4.print_out('Iter.        Energy\n')
        scf_form = '{: >3d}  {: >20.15f}\n'
        i = 0
        converged = False
        e_old = 0
        # Begin SCF iterations
        while not converged and i < self.maxiter:
            F = self.build_fock(D)

            # Transformed Fock matrix to an orthogonal basis
            tF = X*F*X

            # Diagonalized the transformed Fock matrix
            e, tC = la.eigh(tF)

            # Convert eigenvectors to AO basis (from orthogonal basis)
            C = X*tC

            # Form density matrix
            D = C[:, :ndocc]*C[:, :ndocc].T

            # Compute the energy
            E = np.trace((F + T +V)*D) + V_nuc

            i += 1
            # Check convergence
            if abs(E - e_old) < self.e_convergence:
                converged = True
            e_old = E

            psi4.print_out(scf_form.format(i, E))

        psi4.print_out('RHF Energy: {:20.15f}\n'.format(E))

        self.e  = e
        self.C = C
        self.E = E

        return E

    def build_fock(self, D):
        """
        Build J and K - the Coulomb and exchange matrices
        :param D: density matrix
        :return: Fock matrix
        """
        nbf, g = self.nbf, self.g
        J = np.zeros((nbf, nbf))
        K = np.zeros((nbf, nbf))
        for p in range(nbf):
            for q in range(nbf):
                for r in range(nbf):
                    for s in range(nbf):
                        J[p, q] += D[r, s]*g[p, q, r, s]
                        K[p, q] += D[r, s]*g[p, r, q, s]

        F = self.T + self.V + 2*J - K

        return F
