from uhf import UHF
import numpy as np
from geom import molecule
from mp2 import MP2

class uMP2(MP2):

    def __init__(self, mol, basis, scf_conv=6, scf_maxiter=50):
        MP2.__init__(self, mol, basis, scf_conv, scf_maxiter)

    def energy(self):
        uhf = UHF(self.mol,self.basis,convergence=self.scf_conv,maxiter=self.scf_maxiter)
        uhf_E, C, eps, g = uhf.energy(return_integrals=True)
        M = np.einsum('abcd,ai,bj,ck,dl->ijkl',g,C,C,C,C)
        n_occ = uhf.n_occ
        norbitals = uhf.norbitals
        energy = 0
        for i in range(n_occ):
            for j in range(n_occ):
                for a in range(n_occ,norbitals):
                    for b in range(n_occ,norbitals):
                        energy += ((M[i,j,a,b] - M[i,j,b,a])**2) / (eps[i] + eps[j] - eps[a] - eps[b])
        energy = 0.25 * energy
        return energy, uhf_E

if __name__ == '__main__':
    mol = molecule(units='Bohr')
    mol.to_angstrom()
    #uhf = UHF(mol,basis='STO-3G')
    #uhf.energy()
    ump2 = uMP2(mol, basis='STO-3G',scf_maxiter=50)
    print(ump2.energy())