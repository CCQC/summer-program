import numpy as np

from psi4_helper import get_nocc, get_nbf

import sys
sys.path.insert(0,"../../5/aewiens")
#sys.path.insert(0,"../../5/avcopan")

from uhf import UHF

class CIS:

    def __init__(self,mol,mints):

        self.nbf = get_nbf(mints)
        self.nocc = get_nocc(mol)

        uhf = UHF(mol,mints)
        uhf.get_energy()
        self.E0 = uhf.E
        self.e = uhf.e

        self.Gmo = transform_integrals(uhf.g, uhf.C)

        self.norb = 2*self.nbf
        self.nvirtual = self.norb - self.nocc
        self.ndet = self.nvirtual*self.nocc


    def singles_iterator(self):
        """
        return a list of all (i,a) single excitations  xi_i -> xi_a
        """
        return [(i,a) for i in range(self.nocc) for a in range(self.nocc,self.nocc+self.nvirtual)]

    def cis_energies(self):
        """
        Print out a list of CIS excited states and their energies (Eh)
        Return an ordered list of CIS excitation energies
        """
        ##  rename object variables 
        ndet, E0, e, Gmo = self.ndet, self.E0, self.e, self.Gmo

        excite = self.singles_iterator()

        ##  initialize CIS hamiltonian tH
        tH = np.zeros((ndet,ndet))

        ## build CIS Hamiltonian tH
        for P, (i,a) in enumerate(excite):
            for Q, (j,b) in enumerate(excite):
                tH[P,Q] += Gmo[a,j,i,b] + (e[a] - e[i])*(a==b)*(i==j) 
                
        ##  diagonalize tH
        E, C = np.linalg.eigh(tH)

        print( "{:>6s}{:>15s}".format("State","Energy") )
        for i, en in enumerate(E):
            print("{:6d}  {: >16.11f}".format(i,en) )

        self.E = E
        return E

def transform_integrals(g,C):
  return np.einsum('Pp,Pqrs->pqrs', C,
           np.einsum('Qq,PQrs->Pqrs', C,
             np.einsum('Rr,PQRs->PQrs', C,
               np.einsum('Ss,PQRS->PQRs', C, g))))
