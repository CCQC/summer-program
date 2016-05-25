import numpy as np

import sys
sys.path.insert(0,"../../5/aewiens")

from uhf import UHF

class CIS:

    def __init__(self,mol,mints):

        uhf = UHF(mol,mints)
        uhf.compute_energy()

        nbf = uhf.nbf
        self.nocc = uhf.nocc
        self.nvirtual = 2*nbf - self.nocc
        self.ndet = self.nvirtual*self.nocc

        self.E0 = uhf.E
        self.e = uhf.e
        self.C = uhf.C
        self.g = uhf.g

#        self.TF = uhf.TF


    def transform_integrals(self,g,C):
        """
        :param g: 4D array of 2-electron integrals in AO basis
        :param C: 2D array of MO expansion coefficients for AO basis functions
        Return a 4D array (same as g.size) of 2-electron integrals in MO basis
        """
        return np.einsum('Pp,Pqrs->pqrs', C,
                    np.einsum('Qq,PQrs->Pqrs', C,
                        np.einsum('Rr,PQRs->PQrs', C,
                            np.einsum('Ss,PQRS->PQRs', C, g))))


    def get_singles(self):
        """
        return a list of all (i,a) single excitations  xi_i -> xi_a
        """
        return [(i,a) for i in range(self.nocc) for a in range(self.nocc,self.nocc+self.nvirtual)]


    def cis_energies(self):
        """
        Print a list of CIS excited states and their energies (Eh)
        Return an ordered list of CIS excitation energies
        """
        ##  rename object variables 
        ndet, E0, e, g, C = self.ndet, self.E0, self.e, self.g, self.C

#        TF = self.TF

        Gmo = self.transform_integrals(g,C)

        excite = self.get_singles()

        ##  initialize CIS hamiltonian tH
        tH = np.zeros((ndet,ndet))

        ## build CIS Hamiltonian tH
        for P, (i,a) in enumerate(excite):
            for Q, (j,b) in enumerate(excite):
                tH[P,Q] += Gmo[a,j,i,b] + (e[a] - e[i])*(a==b)*(i==j) 
#                tH[P,Q] += Gmo[a,j,i,b] + TF[a,b]*(i==j) - TF[i,j]*(a==b)
                
        ##  diagonalize tH
        E, C = np.linalg.eigh(tH)

        print( "{:>6s}{:>15s}".format("State","Energy") )
        for i, en in enumerate(E):
            print("{:6d}  {: >16.11f}".format(i,en) )

        self.E = E
        return E

