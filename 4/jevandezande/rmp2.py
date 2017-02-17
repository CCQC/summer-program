import psi4
import numpy as np
from mp2 import MP2


class RMP2(MP2):
    """
    RMP2 class
    """

    def __init__(self, rhf, df_basis_name=''):
        super().__init__(rhf, df_basis_name)

    def energy(self):
        """
        Compute the RMP2 energy
        :return: RMP2 energy
        """
        nocc, ntot, gmo, e = self.nocc, self.gmo, self.e, self.ntot

        Ec = 0.0
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nocc, ntot):
                    for b in range(nocc, ntot):
                        Ec += gmo[i, a, j, b]*(2*gmo[i, a, j, b] - gmo[i, b, j, a])/\
                             (e[i] + e[j] - e[a] - e[b])

        self.Ec = Ec
        self.E_mp2 = Ec + self.E_scf

        print('@MP2 correlation energy: {:15.10f}\n'.format(self.Ec))
        print('@Total MP2 energy: {:15.10f}\n'.format(self.E_mp2))

        return self.E_mp2


if __name__ == "__main__":
    import sys
    sys.path.insert(0, '../../3/jevandezande')
    from rhf import RHF
    rhf = RHF('../../3/jevandezande/Options.ini')
    rhf.energy()
    rmp2 = RMP2(rhf)
    rmp2.energy()
