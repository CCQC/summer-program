import numpy as np
from scipy import linalg as la

import sys
sys.path.insert(0,"../../5/aewiens")
from uhf import UHF

from psi4_helper import get_nocc, get_nbf, get_maxiter

class CIS:

    def __init__(self,mol,mints):
        self.nbf = get_nbf(mints)
        self.norb = 2*self.nbf
        self.nocc = get_nocc(mol)
        self.nvirtual = self.norb - self.nocc
        self.maxiter = get_maxiter
        self.uhf = UHF(mol,mints)
        self.E_uhf = self.uhf.get_energy()
        self.e = self.uhf.e
        self.C = self.uhf.C
        self.g = self.uhf.g

    def cis_energy(self):

        ##  rename object variables 
        norb, nocc, nvirtual, e = self.norb, self.nocc, self.nvirtual, self.e

        Gmo = self.transform_integrals()
        Gmo = Gmo - Gmo.swapaxes(2,3)         ### antisymmetrize

        tH = np.zeros((nocc*nvirtual,nocc*nvirtual))
        for i in range(nocc):
            for a in range(nocc,nvirtual):
                for j in range(nocc):
                    for b in range(nocc,nvirtual):
                        tH[j*nvirtual+a,i*nvirtual+b] = (e[a] - e[i])*np.kron(a,b)*np.kron(i,j) + Gmo[a,j,i,b]
        e_cis, L = np.linalg.eigh(tH)
        return e_cis

    def transform_integrals(self):
        g,C = self.g, self.C
        return np.einsum("Pp,Pqrs->pqrs",C,
                    np.einsum("Qq,PQrs->Pqrs",C,
                        np.einsum("Rr,PQRs->PQrs",C,
                            np.einsum("Ss,PQRS->PQRs",C,g))))
