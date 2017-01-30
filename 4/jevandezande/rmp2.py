import psi4
import numpy as np
from itertools import product


class RMP2:
    """
    RMP2 class
    """

    def __init__(self, rhf):
        self.ndocc = rhf.docc
        self.C = rhf.C
        self.e = rhf.e
        self.gao = rhf.g
        self.E_scf= rhf.E_scf

        #self.transform_integrals()
        self.transform_integrals_itertools()
        #self.transform_integrals_noddy()
        #self.transform_integrals_einsum()

    def energy(self):
        """
        Compute the MP2 energy
        :return: MP2 energy
        """
        ndocc, gmo, e = self.ndocc, self.gmo, self.e

        E = 0.0
        for i in range(ndocc):
            for j in range(ndocc):
                for a in range(ndocc, len(e)):
                    for b in range(ndocc, len(e)):
                        E += gmo[i, a, j, b]*(2*gmo[i, a, j, b] - gmo[i, b, j, a])/\
                             (e[i] + e[j] - e[a] - e[b])

        print('@MP2 correlation energy: {:20.15f}\n'.format(E))
        print('@Total MP2 energy: {:20.15f}\n'.format(E + self.E_scf))

        self.emp2 = E

    def transform_integrals_noddy(self):
        """
        Transform the integrals from the AO basis to the MO basis
        """
        C, gao = self.C, self.gao
        nbf = len(self.e)
        self.gmo = np.zeros(gao.shape)
        # Noddy algorithm, O(N^8)
        for p in range(nbf):
            for q in range(nbf):
                for r in range(nbf):
                    for s in range(nbf):
                        for mu in range(nbf):
                            for nu in range(nbf):
                                for rho in range(nbf):
                                    for sigma in range(nbf):
                                        self.gmo[p, q, r, s] += gao[mu, nu, rho, sigma]*\
                                               C[mu, p]*C[nu, q]*C[rho, r]*C[sigma, s]


    def transform_integrals_einsum(self):
        """
        Einsum version of integral transformation
        """
        C, gao = self.C, self.gao
        self.gmo = np.einsum('PQRS,Pp,Qq,Rr,Ss->pqrs', gao, C, C, C, C)

    def transform_integrals(self):
        """
        Transform integrals the efficient way O(N^5)
        """
        C, gao = self.C, self.gao
        nbf = len(self.e)

        g = np.zeros(gao.shape)
        for mu in range(nbf):
            for nu in range(nbf):
                for rho in range(nbf):
                    for s in range(nbf):
                        for sigma in range(nbf):
                            g[mu, nu, rho, s] += gao[mu, nu, rho, sigma]*C[sigma, s]

        gmo = np.zeros(gao.shape)
        for mu in range(nbf):
            for nu in range(nbf):
                for r in range(nbf):
                    for s in range(nbf):
                        for rho in range(nbf):
                            gmo[mu, nu, r, s] += g[mu, nu, rho, s]*C[rho, r]

        g.fill(0)
        for mu in range(nbf):
            for q in range(nbf):
                for r in range(nbf):
                    for s in range(nbf):
                        for nu in range(nbf):
                            g[mu, q, r, s] += gmo[mu, nu, r, s]*C[nu, q]

        gmo.fill(0)
        for p in range(nbf):
            for q in range(nbf):
                for r in range(nbf):
                    for s in range(nbf):
                        for mu in range(nbf):
                            gmo[p, q, r, s] += g[mu, q, r, s]*C[mu, p]

        self.gmo = gmo

    def transform_integrals_itertools(self):
        """
        Simplify iteration using a cartesian product
        """
        C, gao = self.C, self.gao
        nbf = len(self.e)
        g1, g2, g3, g4 = np.zeros(gao.shape), np.zeros(gao.shape), np.zeros(gao.shape), np.zeros(gao.shape)

        for μ, ν, ρ, σ, s in product(range(nbf), repeat=5):
            g1[μ, ν, ρ, s] += gao[μ, ν, ρ, σ]*C[σ, s]
        for μ, ν, ρ, r, s in product(range(nbf), repeat=5):
            g2[μ, ν, r, s] += g1[μ, ν, ρ, s]*C[ρ, r]
        for μ, ν, q, r, s in product(range(nbf), repeat=5):
            g3[μ, q, r, s] += g2[μ, ν, r, s]*C[ν, q]
        for μ, p, q, r, s in product(range(nbf), repeat=5):
            g4[p, q, r, s] += g3[μ, q, r, s]*C[μ, p]

        self.gmo = g4


if __name__ == "__main__":
    import sys
    sys.path.insert(0, '../../3/jevandezande')
    from rhf import RHF
    rhf = RHF('../../3/jevandezande/Options.ini')
    rhf.energy()
    rmp2 = RMP2(rhf)
    rmp2.energy()
