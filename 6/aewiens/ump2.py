
import numpy as np
from scipy import linalg as la

from psi4_helper import get_nbf,get_nocc,get_conv,get_maxiter

import sys
sys.path.insert(0,"../../5/aewiens/")
from uhf import UHF

class UMP2:
    
    def __init__(self,mol,mints):
        self.nocc = get_nocc(mol)
        self.nbf = get_nbf(mints)
        self.norb = 2*self.nbf
        self.conv = get_conv()
        self.maxiter = get_maxiter()
        self.uhf = UHF(mol,mints)
        self.E_uhf = self.uhf.get_energy()
        self.e = self.uhf.e
        self.g = self.uhf.g
        self.C = self.uhf.C

    def get_mp2(self):

        # rename class variables
        nocc, nbf, conv, maxiter,norb = self.nocc, self.nbf, self.conv, self.maxiter,self.norb
        E_uhf, e = self.E_uhf, self.e


        ####  4 different algorithms for integral transformation

        #1 
        Gmo = self.transform_integrals_einsum_faster()
        #2  Gmo = self.transform_integrals_einsum()
        #3  Gmo = self.transform_integrals_einsum_noddy()
        #4  Gmo = self.transform_integrals_noddy()


        Gmo = Gmo - Gmo.swapaxes(2,3)                       ####  Antisymmetrize integrals
        Ecorr = 0.0
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nocc,norb):
                    for b in range(nocc,norb):
                        Ecorr += 0.25 * (Gmo[i,j,a,b])**2 / (e[i] + e[j] - e[a] - e[b])
        return E_uhf + Ecorr


    def transform_integrals_einsum_noddy(self):
        g, C = self.g, self.C
        return np.einsum("PQRS,Pp,Qq,Rr,Ss->pqrs",g,C,C,C,C)

    def transform_integrals_noddy(self):
        g, C = self.g, self.C
        Gmo = np.zeros(g.shape)
        for p in range(self.norb):
            for q in range(self.norb):
                for r in range(self.norb):
                    for s in range(self.norb):
                        for P in range(self.norb):
                            for Q in range(self.norb):
                                for R in range(self.norb):
                                    for S in range(self.norb):
                                        Gmo[p,q,r,s] += C[P,p]*C[Q,q]*C[R,r]*C[S,s] * g[P,Q,R,S]
        return Gmo


    def transform_integrals_einsum_faster(self):
        g,C = self.g, self.C
        return np.einsum("Pp,Pqrs->pqrs",C,
                    np.einsum("Qq,PQrs->Pqrs",C,
                        np.einsum("Rr,PQRs->PQrs",C,
                            np.einsum("Ss,PQRS->PQRs",C,g))))


    def transform_integrals_faster(self):
        g,C,norb = self.g, self.C, self.norb

        Gmo1 = np.zeros(g.shape)
        for P in range(norb):
            for Q in range(norb):
                for R in range(norb):
                    for S in range(norb):
                        for s in range(norb):
                            Gmo1[P,Q,R,s] += C[S,s] * g[P,Q,R,S] 

        Gmo2 = np.zeros(g.shape)
        for P in range(norb):
            for Q in range(norb):
                for R in range(norb):
                    for s in range(norb):
                        for r in range(norb):
                            Gmo2[P,Q,r,s] += C[R,r] * Gmo1[P,Q,R,s] 

        Gmo3 = np.zeros(g.shape)
        for P in range(norb):
            for Q in range(norb):
                for r in range(norb):
                    for s in range(norb):
                        for q in range(norb):
                            Gmo3[P,q,r,s] += C[Q,q] * Gmo2[P,Q,r,s] 

        Gmo = np.zeros(g.shape)
        for P in range(norb):
            for q in range(norb):
                for r in range(norb):
                    for s in range(norb):
                        for p in range(norb):
                            Gmo[p,q,r,s] += C[P,p] * Gmo3[P,q,r,s] 

        return Gmo
