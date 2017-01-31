import abc
import sys
sys.path.insert(0, '../../3/jevandezande')
from integral_transform import transform_integrals_itertools as transform_integrals


class MP2:
    """
    MP2 class
    """

    def __init__(self, scf):
        self.nocc, self.ntot, self.e, self.gao, self.E_scf = scf.nocc, scf.ntot, scf.e, scf.g, scf.E_scf
        self.gmo = transform_integrals(scf.C, scf.g)

    @abc.abstractmethod
    def energy(self):
        """
        Compute the MP2 energy
        :return: MP2 energy
        """

