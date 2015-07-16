import psi4
import numpy as np
import sys
sys.path.append("../../5/jevandezande/")
from uhf import UHF


class MP2:
    """
    Spinorbital version of MP2 (i.e. MP2 using an unrestricted reference)
    """

    def __init__(self, mol, mints):
        self.uhf = UHF(mol, mints)
        self.uhf.compute_energy()
        self.gmo = transform_tei(self.uhf.g, self.uhf.C)

    def compute_energy(self):
        uhf = self.uhf
        e, nocc, norb, gmo = self.uhf.e, self.uhf.nocc, self.uhf.norb, self.gmo

        Ec = 0.0
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nocc, norb):
                    for b in range(nocc, norb):
                        Ec += (1/4.0) * gmo[i, j, a, b]**2 / (e[i]+e[j]-e[a]-e[b])

        self.E = uhf.E + Ec

        psi4.print_out('MP2            E                  Ec\n')
        psi4.print_out('     {:20.15f}  {:20.15f}'.format(self.E, Ec))

        return Ec


def transform_tei(gao, C):
    # g_pqrs = sum_P C_Pp (sum_Q C_Qq (sum_R C_Rr (sum_S C_Ss gao_PQRS)))
    return np.einsum('Pp,Pqrs->pqrs', C,
             np.einsum('Qq,PQrs->Pqrs', C,
               np.einsum('Rr,PQRs->PQrs', C,
                 np.einsum('Ss,PQRS->PQRs', C, gao))))


# Another way of doing the integral transformation:
def transform_tei_einsum(gao, C):
    # g_pqrs = sum_PPQRS gao_PQRS C_Pp C_Qq C_Rr C_Ss
    return np.einsum('PQRS,Pp,Qq,Rr,Ss->pqrs', gao, C, C, C, C)
