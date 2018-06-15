import psi4
import numpy as np
import scipy.linalg as la

class RHF:
    """
    Restricted Hartree-Fock class to be used with PSI4
    """

    def __init__(self, mol, mints):
        """
        Initialize rhf
        :param mol: a psi4 molecule object
        :param mints: a molecular integrals object (from MintsHelper)
        """
        self.mol = mol
        self.mints = mints

        self.V_nuc = mol.nuclear_repulsion_energy()
        self.T = np.matrix(mints.ao_kinetic())
        self.S = np.matrix(mints.ao_overlap())
        self.V = np.matrix(mints.ao_potential())

        self.g = np.array(mints.ao_erl())

        # Determine the number of electrons and the number of doubly occupied orbitals
        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        if mol.multiplicity() != 1 or self.nelec % 2:
            raise Exception("This code only allows closed-shell systems")
        self.ndocc = self.nelec / 2

        self.maxiter = psi4.get_global_option('MAXITER')
        self.e_convergence = psi4.get_global_option('E_CONVERGENCE')

        self.nbf = mins.basisset().nbf()

    def compute_energy(self):
        """
        Compute the rhf energy
        :return: energy
        """
        V_nuc, T, S, V, g = self.V_nuc, self.T, self.S, self.V, self.g
        nbf, ndocc = self.nbf, self.ndocc

        
