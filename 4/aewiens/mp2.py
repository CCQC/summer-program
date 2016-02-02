import numpy as np 
from psi4_helper import get_docc, get_nbf, get_conv, get_maxiter

import sys
sys.path.insert(0,"../../3/aewiens")
from rhf import RHF

class MP2:

    def __init__(self,mol,mints):

        self.docc = get_docc(mol)
        self.nbf = get_nbf(mints)
        self.conv = get_conv()
        self.maxiter = get_maxiter()

        rhf = RHF(mol,mints)
        self.E_scf = rhf.get_energy()

        self.e = rhf.e
        self.G = rhf.G
        self.C = rhf.C


    def get_energy(self):

        #rename class variables
        docc, nbf, conv, maxiter = self.docc, self.nbf, self.conv, self.maxiter
        E_scf, e = self.E_scf, self.e

        Gmo = self.transform_integrals()
        Ecorr = 0.0
        for i in range(docc):
            for j in range(docc):
                for a in range(docc,nbf):
                    for b in range(docc,nbf):
                        Ecorr += (2*Gmo[i,j,a,b]-Gmo[i,j,b,a])*Gmo[i,j,a,b]/(e[i]+e[j]-e[a]-e[b])

        return E_scf + Ecorr

    def transform_integrals(self):
        G,C = self.G, self.C
        return np.einsum("Pqrs,Pp->pqrs",
                np.einsum("PQrs,Qq->Pqrs",
                    np.einsum("PQRs,Rr->PQrs",
                        np.einsum("PQRS,Ss->PQRs", G, C), C) ,C) ,C)


"""
Slower but instructive ways to transform integrals

    def transform_integrals_einsum_2(self):
        G,C = self.G, self.C
        return np.einsum("PQRS,Pp,Qq,Rr,Ss->pqrs",G,C,C,C,C)


    def transform_integrals_noddy(self):
        nbf = self.nbf
        C,G = self.C, self.G

        Gmo = np.zeros(self.G.shape)
        for p in range(nbf):
            for q in range(nbf):
                for r in range(nbf):
                    for s in range(nbf):
                        for i in range(nbf):
                            for j in range(nbf):
                                for k in range(nbf):
                                    for l in range(nbf):
                                        Gmo[p,q,r,s] += C[i,p]*C[j,q]*C[k,r]*C[l,s]*G[i,j,k,l]
        return Gmo
"""
