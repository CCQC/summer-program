import numpy as np
import sys
sys.path.insert(0,"../../../5/aewiens/")
from uhf import UHF

class CEPA0:

    def __init__(self,mol,mints):

        uhf = UHF(mol,mints)
        self.Ehf = uhf.compute_energy()

        self.conv = uhf.conv
        self.maxiter = uhf.maxiter
        
        self.G = uhf.g
        self.C = uhf.C
        self.e = uhf.e
        self.nocc = uhf.nocc
        self.nvirt = uhf.norb - self.nocc

        self.t = np.zeros((self.nocc,self.nocc,self.nvirt,self.nvirt))
        self.Ec = 0.0 


    def transform_integrals(self,g,C):
        return np.einsum("Pp,Pqrs->pqrs",C,
                    np.einsum("Qq,PQrs->Pqrs",C,
                        np.einsum("Rr,PQRs->PQrs",C,
                            np.einsum("Ss,PQRS->PQRs",C,g))))


    def compute_energy(self):
        """
        :Return: CCD correlation energy
        """
        nocc, nvirt, e, t = self.nocc, self.nvirt, self.e, self.t
        g = self.transform_integrals(self.G,self.C)

        o = slice(0,nocc)
        v = slice(nocc,nocc+nvirt)
        x = np.newaxis
        D = e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v] 
        Ep = 1.0/D

        for i in range(self.maxiter):
            ##  terms
            t1 = g[o,o,v,v]
            t2 = np.einsum("abcd,ijcd->ijab",g[v,v,v,v],t)
            t3 = np.einsum("klij,klab->ijab",g[o,o,o,o],t)
            t4 = np.einsum("akic,jkbc->ijab",g[v,o,o,v],t)

            ##  permute t4
            t4P = t4 - t4.transpose((1,0,2,3)) - t4.transpose((0,1,3,2)) + t4.transpose((1,0,3,2))

            ##  update t
            t = Ep*( t1 + 0.5*t2 + 0.5*t3 + t4P )
          
            ##  evaluate Ecorr
            Ec = 0.25 * np.sum(g[o,o,v,v] * t )
            dE = np.fabs(Ec - self.Ec)

            ##  save
            self.Ec = Ec
            self.t = t

            print       ("@CEPA {:3d} {:20.15f} {:20.15f}"  .format(i,Ec,dE) )

            if dE < self.conv: break

        return self.Ec

