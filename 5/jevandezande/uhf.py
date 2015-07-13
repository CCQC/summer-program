import psi4
import numpy as np
import scipy.linalg as la


class UHF:
    """
    Unrestricted Hartree-Fock class for obtaining the unrestricted Hartree-Fock
    energy
    """

    def __init__(self, mol, mints):
        """
        Initialize the uhf
        :param mol: a psi4 molecule object
        :param mints: a molecular integrals object (from MintsHelper)
        """
        self.mol = mol
        self.mints = mints

        self.V_nuc = mol.nuclear_repulsion_energy()
        self.T = block_oei(mints.ao_kinetic())
        self.S = block_oei(mints.ao_overlap())
        self.V = block_oei(mints.ao_potential())

        G = block_tei(mints.ao_eri())
        # Convert to antisymmetrized integrals in physicist's notation
        self.g = G.transpose(0, 2, 1, 3) - G.transpose(0, 2, 3, 1)

        self.h = self.T + self.V

        # S^{-1/2}
        self.X = la.inv(la.sqrtm(self.S)).view(np.matrix)

        # Determine the number of electrons
        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        self.nocc = self.nelec

        self.maxiter = psi4.get_global_option('MAXITER')
        self.e_convergence = psi4.get_global_option('E_CONVERGENCE')

        self.nbf = mints.basisset().nbf()

    def compute_energy(self):
        """
        Compute the rhf energy
        :return: energy
        """
        T, V, h, S, X = self.T, self.V, self.h, self.S, self.X
        V_nuc, nocc = self.V_nuc, self.nocc

        # Use an empty density guess
        D = np.zeros(T.shape)

        psi4.print_out('\nIter.        Energy\n')
        scf_form = '{: >3d}  {: >20.15f}\n'
        i = 0
        converged = False
        e_old = 0
        # Begin SCF iterations
        while not converged and i < self.maxiter:
            # build the Fock with the previous density
            f = self.build_fock(D)

            # Transform the Fock matrix to an orthogonal basis
            tf = X*f*X

            # Diagonalize the transformed Fock matrix
            e, tC = la.eigh(tf)

            # Convert eigenvectors to AO basis (from orthogonal basis)
            C = X*tC

            # Form density matrix
            D = C[:, :nocc]*C[:, :nocc].T

            # Compute the energy
            E = np.trace((h/2 + f/2)*D) + V_nuc

            i += 1
            # Check convergence
            if abs(E - e_old) < self.e_convergence:
                converged = True
            e_old = E

            psi4.print_out(scf_form.format(i, E))

        psi4.print_out('RHF Energy: {:20.15f}\n'.format(E))

        self.E = E

        return E

    def build_fock(self, D):
        """
        Build J and K - the Coulomb and exchange matrices
        :param D: density matrix
        :return: Fock matrix
        """
        v = np.zeros(D.shape)
        for p, q, r, s in np.ndindex(self.g.shape):
            v[p, r] += D[s, q]*self.g[p, q, r, s]

        return self.h + v


# spin blocking functions: transform from spatial orbital {x_mu} basis to
# spin-orbital {x_mu alpha, x_mu beta} basis
# block one-electron integrals
def block_oei(a):
    a = np.matrix(a)
    I2 = np.identity(2)
    return np.matrix(np.kron(I2, a))


# block two-electron integrals
# [must be in chemist's notation, (mu nu|rh si)]
def block_tei(T):
    t = np.array(T)
    n = t.shape[0]
    I2 = np.identity(2)
    T = np.zeros((2*n, 2*n, 2*n, 2*n))
    for p in range(n):
        for q in range(n):
            T[p, q] = np.kron(I2, t[p, q])
    T[n:, n:] = T[:n, :n]
    return T

