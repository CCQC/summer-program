import psi4
import numpy as np
import sys
sys.path.insert(0, '../../5/bz3wp')
from uhf import UHF
import scipy.linalg as la

class MP2:
        
    def __init__(self, uhf):
        self.g_ao = uhf.g
        self.C = uhf.C
        self.e = uhf.e
        self.nocc = uhf.nocc
        self.g_mo = np.zeros_like(self.g_ao) 
        self.nto = uhf.nto

    def transform_integrals(self):
        C, g_ao = self.C, self.g_ao
        self.g_mo = np.einsum('pqrs, pP, qQ, rR, sS-> PQRS',g_ao, C, C, C, C)
         
    def get_energy(self):
        e, g_mo, nocc, nto = self.e, self.g_mo, self.nocc, self.nto
        E = 0.0
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nocc, nto):
                    for b in range(nocc, nto):
                        E += ((1/4)*g_mo[i,j,a,b]**2)/(e[i]+e[j]-e[a]-e[b])  

        print('The second order energy correction is {:20.14f}'.format(E)) 
        return E

if __name__=='__main__':

    uhf = UHF('../../5/bz3wp/Options.ini')
    uhf.get_energy()
    mp2 = MP2(uhf)
    mp2.transform_integrals()
    mp2.get_energy()
