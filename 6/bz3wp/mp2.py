import psi4
import numpy as np
import sys
sys.path.insert(0, '../../5/bz3wp')
from uhf import UHF
import scipy.linalg as la

class MP2:
        
    def __init__(self, uhf):
        self.e = uhf.e
        self.nocc = uhf.nocc
        self.nto = uhf.nto
        self.E0 = uhf.E_SCF
    def transform_integrals(self, C, g_ao):
        return np.einsum('pqrs, pP, qQ, rR, sS-> PQRS',g_ao, C, C, C, C)
         
    def get_energy(self):
        g_mo = self.transform_integrals(uhf.C, uhf.g)
        e, E0, nocc, nto = self.e, self.E0, self.nocc, self.nto
        E = 0.0
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nocc, nto):
                    for b in range(nocc, nto):
                        E += ((1/4)*g_mo[i,j,a,b]**2)/(e[i]+e[j]-e[a]-e[b])  

        print('The second order energy correction is {:20.14f}'.format(E)) 
        print('The total MP2 energy is {:20.14f}'.format(E0+E))
        return E

if __name__=='__main__':

    uhf = UHF('../../5/bz3wp/Options.ini')
    uhf.get_energy()
    mp2 = MP2(uhf)
    mp2.get_energy()
