
import psi4
import numpy as np
import sys
import scipy.linalg as la
sys.path.insert(0, '../../5/bz3wp')
from uhf import UHF

class CCD:
        
    def __init__(self, uhf):                             # getting class variables from UHF
        self.e = uhf.e                                   # UHF orbital energies
        self.E0 = uhf.E_SCF                              # UHF energy
        self.nocc = uhf.nocc                             # number of occupied orbitals
        self.nto = uhf.nto                               # total number of spin orbitals
        self.vir = self.nto - self.nocc                  # number of virtual orbitals
        self.ndet = self.nocc*self.vir
        self.C = uhf.C                                   # MO coefficients
    
    
# testing
if __name__=='__main__':

    uhf = UHF('../../5/bz3wp/Options.ini')
    uhf.get_energy()
