import psi4
import numpy as np

from scf import SCF


class RHF(SCF):
    """
    Restricted Hartree-Fock class
    """
    def __init__(self, options_ini):
        super().__init__(options_ini)
        nocc = self.config['DEFAULT']['docc']
        if nocc.isnumeric():
            self.nocc = int(nocc)
        else:
            self.nocc = nocc.split()
            raise Exception('Occupation arrays currently not supported')
        self.ntot = len(self.H)
        self.nvirt = self.ntot - self.nocc
        
    def energy(self):
        """
        Compute the RHF energy
        :return: the RHF energy
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
        for iteration in range(self.options['SCF_MAX_ITER']):
            # Transform Fock
            tF = A @ F @ A
            # Diagonalize Fock
            e, tC = np.linalg.eigh(tF)
            # Construct new SCF eigenvector
            C = A @ tC
            Cocc = C[:, :nocc]
            # Form new density
            D = Cocc @ Cocc.T
            self.densities.append(D)

            # Construct Fock
            J = np.einsum('pqrs,rs->pq', g, D)
            K = np.einsum('prqs,rs->pq', g, D)
            F = H + 2*J - K
            self.focks.append(F)

            E_scf = np.einsum('pq,pq->', F + H, D) + V_nuc
            energies.append(E_scf)
            ΔE = energies[-1] - energies[-2]
            d_norms.append(np.linalg.norm(self.densities[-2] - self.densities[-1]))
            print('{:3d} {:> 20.14f} {:> 1.5E}  {:>1.5E}'.format(iteration, E_scf, ΔE, d_norms[-1]))

            # Check for convergence
            if abs(ΔE) < 1.0e-10 and d_norms[-1] < 1.0e-10:
                break

            if self.options['DIIS'] and self.options['DIIS_START'] < iteration:
                F = self.extrapolate_diis()

        print('\nEnergy: {:> 15.10f}'.format(E_scf))
        self.C, self.e, self.E_scf = C, e, E_scf
        self.energies, self.d_norms = energies, d_norms
        return E_scf


if __name__ == "__main__":
    water = RHF('Options.ini')
    water.energy()
    #water.plot_convergence()
    #water.plot_density_changes()
