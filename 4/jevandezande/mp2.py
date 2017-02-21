import abc
import sys
import numpy as np
from psi4.core import BasisSet
sys.path.insert(0, '../../3/jevandezande')
from integral_transform import transform_integrals_einsum as transform_integrals, transform_integrals_df


class MP2:
    """
    MP2 class
    """

    def __init__(self, scf, df_basis_name=''):
        self.df_basis_name = df_basis_name
        self.nocc, self.ntot, self.e, self.gao, self.E_scf, self.molecule = scf.nocc, scf.ntot, scf.e, scf.g, scf.E_scf, scf.molecule

        if self.df_basis_name:
            df_basis = BasisSet.build(self.molecule, 'DF_BASIS_MP2', df_basis_name, puream=0)
            self.gmo = transform_integrals_df(scf.C, scf.basis, df_basis, bool(scf.spin))
        else:
            if scf.spin != 0 and len(scf.H) == len(self.gao):
                self.gao = block_tei(self.gao)
            self.gmo = transform_integrals(scf.C, self.gao)

    @abc.abstractmethod
    def energy(self):
        """
        Compute the MP2 energy
        :return: MP2 energy
        """
        pass


# Spin blocking functions:
# transform from spatial orbital {x_μ} basis to spin-orbital {x_μα, x_μβ} basis
def block_tei(a):
    """
    Spin block two-electron integrals
    """
    I = np.eye(2)
    a = np.kron(I, a)
    return np.kron(I, a.T)
