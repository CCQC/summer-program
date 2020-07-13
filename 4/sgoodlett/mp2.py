from scf import SCF
from geom import molecule
import psi4
import numpy as np
from scipy.linalg import inv
class MP2():

    def __init__(self,mol,basis,scf_conv=6,scf_max_iter=50):
        self.mol = mol
        self.basis =basis
        self.scf_conv = scf_conv
        self.scf_max_iter = scf_max_iter

    def mp2_energy(self):
        #Run RHF to obtain MO coefficients C and orbital energies
        scf = SCF(self.mol,self.basis,convergence=self.scf_conv,maxiter=self.scf_max_iter)
        scf_E, C, eps, g = scf.energy(return_integrals=True)
        #Transform 2-e ints from AO to MO basis (Eq. 1)
        M = np.einsum('ijkl,ip,jq,kr,ls->pqrs',g,C,C,C,C)
        n_occ = scf.n_occ
        n_virt = scf.norbitals - scf.n_occ
        n_orbitals = scf.norbitals
        #Eval. MP2 energy (Eq. 2)
        mp2_energy = 0
        mp2_tE = 0
        print(eps)
        for i in range(n_occ):
            for j in range(n_occ):
                for a in range(n_occ,n_orbitals):
                    for b in range(n_occ,n_orbitals):
                        mp2_energy += (M[i,j,a,b]*(2*M[i,j,a,b] - M[i,j,b,a])) / (eps[i] + eps[j] - eps[a] - eps[b])
        return mp2_energy, scf_E

if __name__ == '__main__':
    mol = molecule(units='Bohr')
    mol.to_angstrom()
    mp2 = MP2(mol,'STO-3G',scf_conv=9,scf_max_iter=500)
    mp2_E, scf_E = mp2.mp2_energy()
    print('HF = {}, MP2 = {}, Total_E = {}'.format(scf_E,mp2_E,scf_E+mp2_E))