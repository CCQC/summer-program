from __future__ import print_function
import psi4
import numpy as np


class MP2:
    """
    MP2 class
    """

    def __init__(self, rhf_obj):
        self.ndocc = rhf_obj.ndocc
        self.nbf = rhf_obj.nbf
        self.C = rhf_obj.C
        self.e = rhf_obj.e
        self.gao = rhf_obj.g
        self.ehf = rhf_obj.E

        self.transform_integrals()
        #self.transform_integrals_noddy()
        #self.transform_integrals_einsum()

    def compute_energy(self):
        """
        Compute the MP2 energy
        :return: MP2 energy
        """
        ndocc = self.ndocc
        nbf = self.nbf
        gmo = self.gmo
        e = self.e
        E = 0.0
        for i in range(ndocc):
            for j in range(ndocc):
                for a in range(ndocc, nbf):
                    for b in range(ndocc, nbf):
                        E += gmo[i, a, j, b]*(2*gmo[i, a, j, b] - gmo[i, b, j, a])/\
                             (e[i] + e[j] - e[a] - e[b])

        psi4.print_out('@MP2 correlation energy: {:20.15f}\n'.format(E))
        psi4.print_out('@Total MP2 energy: {:20.15f}\n'.format(E + self.ehf))

        self.emp2 = E

    def transform_integrals_noddy(self):
        """
        Transform the integrals from the AO basis to the MO basis
        """
        nbf = self.nbf
        C = self.C
        gao = self.gao
        gmo = np.zeros(gao.shape)
        # Noddy algorithm, O(N^8)
        print(nbf)
        for p in range(nbf):
            print(p, sep=", ")
            for q in range(nbf):
                for r in range(nbf):
                    for s in range(nbf):
                        for mu in range(nbf):
                            for nu in range(nbf):
                                for rho in range(nbf):
                                    for sigma in range(nbf):
                                        gmo[p, q, r, s] += gao[mu, nu, rho, sigma]*\
                                               C[mu, p]*C[nu, q]*C[rho, r]*C[sigma, s]
        print()

        self.gmo = gmo

    def transform_integrals_einsum(self):
        """
        Einsum version of integral transformation
        """
        gao = self.gao
        C = self.C
        self.gmo = np.einsum('PQRS,Pp,Qq,Rr,Ss->pqrs', gao, C, C, C, C)

    def transform_integrals(self):
        """
        Transform integrals the efficient way O(N^5)
        """
        gao = self.gao
        nbf = self.nbf
        C = self.C

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
