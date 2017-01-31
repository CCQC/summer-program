# Configuration Interaction Singles class

import psi4
import numpy as np
import sys
import scipy.linalg as la
sys.path.insert(0, '../../5/bz3wp')
from uhf import UHF

class CIS:
        
    def __init__(self, uhf):                             # getting class variables from UHF
        self.e = uhf.e                                   # UHF orbital energies
        self.E0 = uhf.E_SCF                              # UHF energy
        self.nocc = uhf.nocc                             # number of occupied orbitals
        self.nto = uhf.nto                               # total number of spin orbitals
        self.vir = self.nto - self.nocc                  # number of virtual orbitals
        self.ndet = self.nocc*self.vir
        self.C = uhf.C                                   # MO coefficients
    
    def transform_integrals(self, C, g_ao):              # integral transformation function
        return np.einsum('pqrs, pP, qQ, rR, sS-> PQRS',g_ao, C, C, C, C)
         
    def single_excite(self):                             # all possible single excitations
        return[(i,a) for i in range(self.nocc) for a in range(self.nocc, self.nto)]
 
    def get_energy(self):
        
        ndet,E0, e, C = self.ndet, self.E0, self.e, self.C
        
        g_mo = self.transform_integrals(C, uhf.g)        # call integral transformation function
        
        Hc = np.zeros((ndet, ndet))                      # Hcorrelation = Hcore - E0

        for P, (i,a) in enumerate(self.single_excite()): # calculating matrix elements of single excitations
            for Q, (j,b) in enumerate(self.single_excite()):
                Hc[P,Q] = g_mo[a,j,i,b] + (e[a]-e[i])*(a==b)*(i==j)

        E, C = np.linalg.eigh(Hc)                        # diagonalize correlation hamiltonian

        
        print('{:>6s}{:>10s}'.format('State', 'Energy'))
        for i, ei in enumerate(E): 
            print('{:6d} {:>16.11f}'.format(i,ei))
        
        self.E = E                                       # saving value
        
        return E
    
# testing
if __name__=='__main__':

    uhf = UHF('../../5/bz3wp/Options.ini')
    uhf.get_energy()
    cis = CIS(uhf)
    cis.get_energy()
