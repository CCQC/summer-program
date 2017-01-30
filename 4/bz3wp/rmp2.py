import psi4
import numpy as np
import sys
sys.path.insert(0, '../../3/bz3wp')
from rhf import RHF


class RMP2:
        
    def __init__(self, rhf):
        self.G = rhf.G
        self.C = rhf.C
        self.e = rhf.e
        self.ndocc = rhf.ndocc
        self.ntot = rhf.ntot
        self.E0 = rhf.E_SCF

    def transform_integrals(self, C, G):
        
        return np.einsum('pqrs, pP, qQ, rR, sS-> PQRS', G, C, C, C, C)
         
    def get_energy(self):
        e, E0, ndocc, ntot = self.e, self.E0, self.ndocc, self.ntot
        g_mo = self.transform_integrals(rhf.C, rhf.G)
        E = 0.0
        for i in range(ndocc):
            for j in range(ndocc):
                for a in range(ndocc, ntot):
                    for b in range(ndocc, ntot):        
                        E += (g_mo[i,a,j,b]*(2*g_mo[i,a,j,b]-g_mo[i,b,j,a]))/(e[i]+e[j]-e[a]-e[b])  

        print('The second order energy correction is {:20.14f}'.format(E)) 
        print('The total RMP2 energy is {:20.14f}'.format(E0 + E))
        return E


if __name__=='__main__':

    rhf = RHF('../../3/bz3wp/Options.ini')
    rhf.get_energy()
    rmp2 = RMP2(rhf)
    rmp2.get_energy()
