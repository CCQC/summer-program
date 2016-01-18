import numpy as np

import sys
sys.path.insert(0,"../../5/aewiens")

from uhf import UHF

class CCD:

    def __init__(self,mol,mints):

        uhf = UHF(mol, mints)
        uhf.get_energy()

        # import from uhf
        self.e, self.C, self.g = uhf.e, uhf.C, uhf.g
        self.maxiter = uhf.maxiter
        self.conv = uhf.conv
        nbf = uhf.nbf
        self.nocc = uhf.nocc
        self.nvirt = 2*nbf - self.nocc

        self.Ec = 0.0 
        self.t = np.zeros((self.nvirt,self.nvirt,self.nocc,self.nocc))


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


    def amplitude(self):

        Gmo = self.transform_integrals(self.g, self.C)
        t, nocc = self.t, self.nocc

        return Gmo[nocc:,nocc:,:nocc,:nocc] + 0.5*np.einsum("abcd,cdij->abij",Gmo[nocc:,nocc:,nocc:,nocc:],t) + 0.5*np.einsum("klij,abkl->abij",Gmo[:nocc,:nocc,:nocc,:nocc],t) 
        #+ np.einsum("acik,bcjk->abij",Gmo[nocc:,:nocc,:nocc,nocc:],t) 
        #- np.einsum("akjc,bkic->aijb",Gmo[nocc:,:nocc,:nocc,nocc:],t) - np.einsum("bkic,akjc->bjia",Gmo[nocc:,:nocc,:nocc,nocc:],t) - np.einsum("bkjc,akic->bija",Gmo[nocc:,:nocc,:nocc,nocc:],t)


    def win(self):

        for k in range(self.maxiter):
            Gmo = self.transform_integrals(self.g,self.C)

            t, nocc = self.t, self.nocc

            t_old, Ec_old = self.t, self.Ec
            
            t_new = np.multiply(t_old,self.amplitude())
            
            Ec_new = np.einsum("cdkl,cdkl",Gmo[nocc:,nocc:,:nocc,:nocc],t)

            dE = np.fabs(Ec_new - Ec_old)

            # save stuff
            self.Ec = Ec_new
            self.t = t_new

            if dE < self.conv: break

        return self.Ec

