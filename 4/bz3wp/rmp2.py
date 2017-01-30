import psi4
import numpy as np
import sys
sys.path.insert(0, '../../3/bz3wp')
from rhf import RHF


class RMP2:
        
    def __init__(self, rhf):
        self.I = rhf.I
        self.C = rhf.C
        self.e = rhf.e
        self.ndocc = rhf.ndocc
        self.ntot = rhf.ntot
        self.g = np.zeros_like(self.I) 

    def transform_integrals(self):
        
        C, I = self.C, self.I
        self.g = np.einsum('pqrs, pP, qQ, rR, sS-> PQRS',I, C, C, C, C)
         
    def get_energy(self):
        e, g, ndocc, ntot = self.e, self.g, self.ndocc, self.ntot
        E = 0.0
        for i in range(ndocc):
            for j in range(ndocc):
                for a in range(ndocc, ntot):
                    for b in range(ndocc, ntot):        
                        E += (g[i,a,j,b]*(2*g[i,a,j,b]-g[i,b,j,a]))/(e[i]+e[j]-e[a]-e[b])  

        print('The second order energy correction is {:20.14f}'.format(E)) 
        return E

if __name__=='__main__':

    rhf = RHF('../../3/bz3wp/Options.ini')
    rhf.get_energy()
    rmp2 = RMP2(rhf)
    rmp2.transform_integrals()
    rmp2.get_energy()
