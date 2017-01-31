# Open-shell MP2 with UHF reference class

import psi4
import numpy as np
import sys
sys.path.insert(0, '../../5/bz3wp')
from uhf import UHF
import scipy.linalg as la

class MP2:
        
    def __init__(self, uhf):
        self.e = uhf.e                      # gets orbital energies from UHF
        self.nocc = uhf.nocc                # gets number off occupied orbitals from UHF
        self.nto = uhf.nto                  # gets total number of spin orbitals from UHF
        self.E0 = uhf.E_SCF                 # UHF energy

    def transform_integrals(self, C, g_ao): #Transforms from AO to MO
        return np.einsum('pqrs, pP, qQ, rR, sS-> PQRS',g_ao, C, C, C, C)
         
    def get_energy(self):
        # calls transform_integrals function on uhf coeffiecients and two-electron integral
        g_mo = self.transform_integrals(uhf.C, uhf.g)
        e, E0, nocc, nto = self.e, self.E0, self.nocc, self.nto
        E = 0.0
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nocc, nto):
                    for b in range(nocc, nto):
                        E += ((1/4)*g_mo[i,j,a,b]**2)/(e[i]+e[j]-e[a]-e[b])  # MP2 energy expression

        print('The second order energy correction is {:20.14f}'.format(E)) 
        print('The total MP2 energy is {:20.14f}'.format(E0+E))
        return E

# testing
if __name__=='__main__':

    uhf = UHF('../../5/bz3wp/Options.ini')
    uhf.get_energy()
    mp2 = MP2(uhf)
    mp2.get_energy()
