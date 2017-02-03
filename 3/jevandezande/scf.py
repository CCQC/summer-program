import psi4
import numpy as np
import configparser
import numpy as np
from matplotlib import pyplot as plt
import abc


class SCF:
    """
    A stub SCF class to be used for RHF, UHF, etc.
    """
    def __init__(self, options_ini):
        self.config = configparser.ConfigParser()
        self.config.read(options_ini)
        self.molecule = psi4.geometry(self.config['DEFAULT']['molecule'])
        self.molecule.update_geometry()
        self.V_nuc = self.molecule.nuclear_repulsion_energy()

        self.options = {}
        self.options['BASIS'] = self.config['DEFAULT']['basis']
        self.options['SCF_MAX_ITER'] = self.config.getint('SCF', 'max_iter', fallback=50)

        self.options['DIIS'] = self.config.getboolean('SCF', 'diis', fallback=True)
        self.options['DIIS_NVECTOR'] = self.config.getint('SCF', 'diis_nvector', fallback=6)
        self.options['DIIS_START'] = self.config.getint('SCF', 'diis_start', fallback=6)

        self.basis = psi4.core.BasisSet.build(self.molecule, 'BASIS', self.options['BASIS'], puream=0)
        mints = psi4.core.MintsHelper(self.basis)

        self.S = mints.ao_overlap().to_array()
        self.T = mints.ao_kinetic().to_array()
        self.V = mints.ao_potential().to_array()
        self.g = mints.ao_eri().to_array()

        self.H = self.T + self.V

        A = mints.ao_overlap()
        A.power(-0.5, 1.e-16)
        self.A = A.to_array()

    @abc.abstractmethod
    def energy(self):
        """
        Compute the SCF energy
        """
        return

    def extrapolate_diis(self):
        """
        Extrapolate the Fock matrix
        e = FDS - SDF
        P q = f
        +---+---+---+---+ +---+   +---+
        |P00|P01|P02|-1 | | q0|   | 0 |
        +---+---+---+---+ +---+   +---+
        |P10|P11|P12|-1 | | q1|   | 0 |
        +---+---+---+---+ +---+ = +---+
        |P20|P21|P22|-1 | | q2|   | 0 |
        +---+---+---+---+ +---+   +---+
        |-1 |-1 |-1 | 0 | | λ |   |-1 |
        +---+---+---+---+ +---+   +---+
        """
        e_vec = []
        S = self.S
        start = max(1, len(self.focks) - self.options['DIIS_NVECTOR'])
        num_mats = len(self.focks[start:])
        for F, D in zip(self.focks[start:], self.densities[start:]):
            e_vec.append(F @ D @ S - S @ D @ F)
        
        P = np.zeros((num_mats + 1, num_mats + 1))
        f = np.zeros((num_mats + 1))
        f[-1] = -1

        for i in range(num_mats):
            P[-1,  i] = -1
            P[ i, -1] = -1
            P[i, i] = (e_vec[i]**2).sum()
            for j in range(i + 1, num_mats):
                P[i, j] = P[j, i] = (e_vec[i]*e_vec[j]).sum()

        q_vec = np.linalg.solve(P, f)

        F_diis = np.zeros_like(self.focks[0])
        for q, F in zip(q_vec[:-1], self.focks[start:]):
            F_diis += q*F

        return F_diis

    def plot_convergence(self):
        """
        Plot the convergence of energy change and density norm
        """
        energies, d_norms = np.array(self.energies), self.d_norms
        e, = plt.plot(abs(energies[1:-1] - energies[2:]), label=r'$\Delta$ E')
        d, = plt.plot(d_norms[1:], label=r'$||\Delta$ D$||$')
        plt.yscale('log')
        plt.legend([e, d])
        plt.show()

    def plot_densities(self):
        """
        Plot the densities
        """
        densities = self.densities
        fig, axes = plt.subplots(3, 3)
        axes = axes.reshape(-1)
        cmap = plt.get_cmap('Oranges_r')

        for i in range(1, min(len(densities), 10)):
            # Select distributed throughout
            j = i
            d = densities[j]
            if len(densities) > 10:
                j = len(densities)//10 * i
                d = densities[j]
            ax = axes[i - 1]
            ax.set_title('{:d}'.format(j))
            ax.imshow(d, interpolation='nearest', cmap=cmap)

        #fig.colorbar(densities[-1])
        plt.show()


    def plot_density_changes(self):
        """
        Plot the change in density
        """
        densities = self.densities
        fig, axes = plt.subplots(3, 3)
        axes = axes.reshape(-1)
        cmap = plt.get_cmap('Oranges_r')

        ims = []
        for i in range(1, min(len(densities), 10)):
            # Select distributed throughout
            j = i
            d_change = densities[j] - densities[j-1]
            if len(densities) > 10:
                j = len(densities)//10 * i
                d = densities[j] - densities[j-1]
            ax = axes[i - 1]
            ax.set_title('{:d}'.format(j))
            ims.append(ax.imshow(d, interpolation='nearest', cmap=cmap))

        fig.colorbar(ims[-1])
        plt.show()
