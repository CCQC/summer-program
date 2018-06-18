import psi4
import numpy as np
import abc
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../../3/jevandezande')
from integral_transform import transform_integrals_itertools as transform_integrals


class CC:
    
    def __init__(self, uhf):
        self.uhf = uhf
        self.ntot, self.nocc, self.nvirt = uhf.ntot, uhf.nocc, uhf.nvirt
        self.gmo = transform_integrals(uhf.C, uhf.g)
        self.maxiter = 50
        self.E_convergence = 1e-5
        self.method_name = 'CC'

    def energy(self):
        """
        Compute the coupled cluster energy
        """
        energies = [0.0]
        nocc, nvirt, g = self.nocc, self.nvirt, self.gmo
        t  = np.zeros((nocc, nocc, nvirt, nvirt))
        o  = slice(None, nocc)
        v  = slice(nocc, None)
        e = self.uhf.e
        x = np.newaxis
        Ep = 1./(e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v])

        print('{}'.format(self.method_name))
        print('Iter.     Energy       ΔE')
        print('-----------------------------')
        for i in range(self.maxiter):
            t = self.update_amplitudes(t)
            t *= Ep
            E_cc = 1/4*np.sum(g[o, o, v, v] * t)
            # evaluate energy
            E_cc = 1./4 * np.sum(g[o,o,v,v] * t)
            energies.append(E_cc)
            ΔE = energies[-1] - energies[-2]
            print('{: >3d} {: >15.10f} {: >12.5E}'.format(i, E_cc, ΔE))
            if abs(ΔE) < self.E_convergence:
                break

        self.E_cc, self.t, self.energies = E_cc, t, energies
        return E_cc

    @abc.abstractmethod
    def update_amplitudes(self, t):
        """
        Update the CC amplitudes (differs for each theory level)
        :return: new t amplitudes
        """

    def plot_convergence(self):
        """
        Plot the convergence of energy change and density norm
        """
        energies = np.array(self.energies)
        e, = plt.plot(abs(energies[1:-1] - energies[2:]), label=r'$\Delta$ E')
        plt.yscale('log')
        plt.legend([e])
        plt.show()


def permutations(w, *args):
    """
    Generates all permutations of the arguments
    """
