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

    def compute_energy(self):
        """
        Compute the rhf energy
        :return: energy
        """
        pass

