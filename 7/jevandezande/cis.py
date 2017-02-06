import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,"../../3/jevandezande")
from integral_transform import transform_integrals_itertools as transform_integrals


class CIS:

    def __init__(self, uhf):
        self.uhf = uhf
        self.nocc = uhf.nocc
        self.ntot = len(uhf.H)
        self.nvirt = self.ntot - self.nocc
        self.ndet = self.nocc*self.nvirt
        self.gmo = transform_integrals(uhf.C, uhf.g)

    def get_singles(self):
        """
        Generate all (i,a) single excitations  i -> a
        """
        for i in range(self.nocc):
            for a in range(self.nocc, self.ntot):
                yield (i, a)

    def energy(self):
        """
        CIS energy
        :return: Lowest CIS energy (all energies are stored in self.energies_cis)
        """
        ndet, e, gmo = self.ndet, self.uhf.e, self.gmo

        #  Make empty Hamiltonian
        H = np.zeros((ndet, ndet))

        # Build Hamiltonian
        for P, (i, a) in enumerate(self.get_singles()):
            for Q, (j, b) in enumerate(self.get_singles()):
                # Slater's rules: h_ia + g_ijab
                H[P, Q] = (e[a] - e[i])*(a==b)*(i==j) + gmo[a, j, i, b]

        e_vec, C = np.linalg.eigh(H)

        print('State      Energy             ΔE')
        print('-------------------------------------')
        for i, en in enumerate(e_vec[:10]):
            print('{: >5d} {: >15.10f} {: >15.10f}'.format(i, self.uhf.E_scf + en, en))
        if len(e_vec) > 10:
            print('...\n')

        self.e_cis, self.energies_cis, self.C, self.H = e_vec[0], e_vec, C, H
        return self.e_cis

    def show_hamiltonian(self):
        """
        Display the hamiltonian
        """
        H = self.H
        cmap = plt.get_cmap('viridis')
        im = plt.matshow(H, interpolation='nearest', cmap=cmap, norm=LogNorm(vmin=1.0e-10, vmax=H.max()))
        plt.colorbar(im)
        plt.show()


if __name__ == "__main__":
    sys.path.insert(0, '../../5/jevandezande')
    from uhf import UHF
    uhf = UHF('../../3/jevandezande/Options.ini')
    uhf.energy()
    cis = CIS(uhf)
    cis.energy()
    cis.show_hamiltonian()
