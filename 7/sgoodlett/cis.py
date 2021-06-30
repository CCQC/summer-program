from uhf import UHF
from geom import molecule
import numpy as np
from itertools import product
from scipy.linalg import eigh

def kron_delta(a,b):
    if a == b:
        return 1
    else:
        return 0

class CIS():

    def __init__(self,uhf):
        self.uhf = uhf
        self.uhf_E, self.C, self.eps, self.g = self.uhf.energy(return_integrals=True)

    def energy(self):
        M = np.einsum('abcd,ap,bq,cr,ds->pqrs',self.g - self.g.transpose((0,1,3,2)),self.C,self.C,self.C,self.C)
        singles = []
        for i in range(self.uhf.n_occ):
            for a in range(self.uhf.n_occ,self.uhf.norbitals):
                singles.append((i,a))
        H = np.zeros((len(singles),len(singles)))
        for p in range(len(singles)):
            for q in range(len(singles)):
                i,a = singles[p]
                j,b = singles[q]
                H[p][q] = (self.eps[a] - self.eps[i]) * kron_delta(i,j) * kron_delta(a,b) + M[a,j,i,b]
        E_K, c_K = eigh(H)
        return E_K


if __name__ == '__main__':
    mol = molecule(units='Bohr')
    mol.to_angstrom()
    uhf = UHF(mol,basis='STO-3G')
    cis = CIS(uhf)
    print(cis.energy())