import psi4
import numpy as np
import scipy.linalg as la
import sys
sys.path.insert(0, '../../3/jevandezande/')
from scf import SCF


class UHF(SCF):
    """
    Unrestricted Hartree-Fock class
    """

    def __init__(self, options_ini):
        super().__init__(options_ini)

        self.S = block_oei(self.S)
        self.T = block_oei(self.T)
        self.V = block_oei(self.V)
        # Convert to antisymmetrized integrals in physicist's notation
        G = block_tei(self.g)
        self.g = G.transpose(0, 2, 1, 3) - G.transpose(0, 2, 3, 1)

        self.H = self.T + self.V

        # S^{-1/2}
        self.A = la.inv(la.sqrtm(self.S))

        # Determine the number of electrons
        mol = self.molecule
        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        self.ntot = len(self.H)
        self.nocc = self.nelec
        self.nvirt = self.ntot - self.nocc

        self.norb = len(self.H)

    def energy(self):
        """
        Compute the uhf energy
        :return: energy
        """
        energies = [0.0]
        self.densities = [np.zeros_like(self.H)]
        self.focks = [np.zeros_like(self.H)]
        d_norms = []

        g, H, A, nocc, V_nuc = self.g, self.H, self.A, self.nocc, self.V_nuc

        # Core guess
        F = H

        print('Iter         Energy            ΔE         ‖ΔD‖')
        print('--------------------------------------------------')
        # Begin SCF iterations
        for iteration in range(self.options['SCF_MAX_ITER']):
            # Transform Fock matrix
            tF = A @ F @ A

            # Diagonalize Fock matrix
            e, tC = la.eigh(tF)

            # Construct new SCF eigenvectors
            C = A @ tC

            # Form new density
            D = C[:, :nocc] @ C[:, :nocc].T
            self.densities.append(D)

            # Construct Fock
            F = self.build_fock(D)
            self.focks.append(F)

            # Compute energy
            E_scf = np.trace((H + F) @ D)/2 + V_nuc

            energies.append(E_scf)
            ΔE = energies[-1] - energies[-2]
            d_norms.append(np.linalg.norm(self.densities[-2] - self.densities[-1]))
            print('{:3d} {:> 20.14f} {:> 1.5E}  {:>1.5E}'.format(iteration, E_scf, ΔE, d_norms[-1]))

            if abs(ΔE) < 1.0e-10 and d_norms[-1] < 1.0e-10:
                break

            if self.options['DIIS'] and self.options['DIIS_START'] < iteration:
                F = self.extrapolate_diis()

        print('\nUHF Energy: {:15.10f}\n'.format(E_scf))
        self.energies, self.d_norms = energies, d_norms
        self.C, self.e, self.E_scf = C, e, E_scf

        return E_scf

    def build_fock(self, D):
        """
        Build J and K - the Coulomb and exchange matrices
        :param D: density matrix
        :return: Fock matrix
        """
        v = np.zeros_like(D)
        for p, q, r, s in np.ndindex(self.g.shape):
            v[p, r] += D[s, q]*self.g[p, q, r, s]

        return self.H + v


# Spin blocking functions:
# transform from spatial orbital {x_μ} basis to spin-orbital {x_μα, x_μβ} basis
def block_oei(a):
    """
    Spin block one-electron integrals
    """
    return la.block_diag(a, a)


def block_tei(a):
    """
    Spin block two-electron integrals
    """
    I = np.eye(2)
    a = np.kron(I, a)
    return np.kron(I, a.T)


if __name__ == "__main__":
    water = UHF('../../3/jevandezande/Options.ini')
    water.energy()
    water.plot_convergence()
    water.plot_density_changes()
