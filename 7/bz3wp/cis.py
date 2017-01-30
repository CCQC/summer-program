import psi4
import numpy as np
import sys
sys.path.insert(0, '../../5/bz3wp')
from uhf import UHF
import scipy.linalg as la

class CIS:
        
    def __init__(self, uhf):
        self.e = uhf.e
        self.E0 = uhf.E_SCF
        self.nocc = uhf.nocc
        self.nto = uhf.nto
        self.vir = self.nto - self.nocc
        self.ndet = self.nocc*self.vir
        self.C = uhf.C
    def transform_integrals(self, C, g_ao):
        return np.einsum('pqrs, pP, qQ, rR, sS-> PQRS',g_ao, C, C, C, C)
         
    def single_excite(self):
        return[(i,a) for i in range(self.nocc) for a in range(self.nocc, self.nto)]
 
    def get_energy(self):
        ndet,E0, e, C = self.ndet, self.E0, self.e, self.C
        g_mo = self.transform_integrals(C, uhf.g)
        
        Hc = np.zeros((ndet, ndet))

        for P, (i,a) in enumerate(self.single_excite()):
            for Q, (j,b) in enumerate(self.single_excite()):
                Hc[P,Q] = g_mo[a,j,i,b] + (e[a]-e[i])*(a==b)*(i==j)
        E, C = np.linalg.eigh(Hc)

        
        print('{:>6s}{:>10s}'.format('State', 'Energy'))
        for i, ei in enumerate(E): 
            print('{:6d} {:>16.11f}'.format(i,ei))
        self.E =E    
        return E
    
if __name__=='__main__':

    uhf = UHF('../../5/bz3wp/Options.ini')
    uhf.get_energy()
    cis = CIS(uhf)
    cis.get_energy()
