import numpy as np
import sys
sys.path.insert(0,"../../5/aewiens")

from uhf import UHF


class CCD:

    def __init__(self,mol,mints):

        uhf = UHF(mol, mints)
        uhf.compute_energy()

        # import from uhf
        self.e, self.C, self.G = uhf.e, uhf.C, uhf.g
        self.maxiter, self.conv = uhf.maxiter, uhf.conv

        self.nbf, self.nocc = uhf.nbf, uhf.nocc
        self.nvirt = 2*self.nbf - self.nocc

        self.Ec = 0.0 
        self.t = np.zeros((self.nocc,self.nocc,self.nvirt,self.nvirt))


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

    def get_amplitudes(self):

        g = self.transform_integrals(self.G, self.C)
        nocc, t, e, = self.nocc, self.t, self.e

        o = slice(None,nocc)
        v = slice(nocc,None)
        x = np.newaxis
        Ep = 1.0/ (e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v] )

        for k in range(self.maxiter):
            # terms
            t1 = g[o,o,v,v]
            t2 = np.einsum("abcd,ijcd->ijab",g[v,v,v,v],t)
            t3 = np.einsum("klij,klab->ijab",g[o,o,o,o],t)
            t4 = np.einsum("akic,jkbc->ijab",g[v,o,o,v],t)
            t5 = np.einsum("klcd,ijac,klbd->ijab",g[o,o,v,v],t,t)
            t6 = np.einsum("klcd,ikab,jlcd->ijab",g[o,o,v,v],t,t)
            t7 = np.einsum("klcd,ijcd,klab->ijab",g[o,o,v,v],t,t)
            t8 = np.einsum("klcd,ikac,jlbd->ijab",g[o,o,v,v],t,t)

            # permutations
            t4P = t4 - t4.transpose((1,0,2,3)) - t4.transpose((0,1,3,2)) + t4.transpose((1,0,3,2))
            t5P = t5 - t5.transpose((0,1,3,2))
            t6P = t6 - t6.transpose((1,0,2,3))
            t8P = t8 - t8.transpose((1,0,2,3))

            # update t2 amplitudes
            t = t1 + 0.5*t2 + 0.5*t3 + t4P - 0.5*t5P - 0.5*t6P + 0.25*t7 + t8P
            t *= Ep 

        return t
    
    def compute_energy(self):


    def get_energy(self):
        """
        :Return: CCD correlation energy
        """

        g = self.transform_integrals(self.G, self.C)
        nocc, t, e, = self.nocc, self.t, self.e

        o = slice(None,nocc)
        v = slice(nocc,None)
        x = np.newaxis
        Ep = 1.0/ (e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v] )

        for k in range(self.maxiter):
            # terms
            t1 = g[o,o,v,v]
            t2 = np.einsum("abcd,ijcd->ijab",g[v,v,v,v],t)
            t3 = np.einsum("klij,klab->ijab",g[o,o,o,o],t)
            t4 = np.einsum("akic,jkbc->ijab",g[v,o,o,v],t)
            t5 = np.einsum("klcd,ijac,klbd->ijab",g[o,o,v,v],t,t)
            t6 = np.einsum("klcd,ikab,jlcd->ijab",g[o,o,v,v],t,t)
            t7 = np.einsum("klcd,ijcd,klab->ijab",g[o,o,v,v],t,t)
            t8 = np.einsum("klcd,ikac,jlbd->ijab",g[o,o,v,v],t,t)

            # permutations
            t4P = t4 - t4.transpose((1,0,2,3)) - t4.transpose((0,1,3,2)) + t4.transpose((1,0,3,2))
            t5P = t5 - t5.transpose((0,1,3,2))
            t6P = t6 - t6.transpose((1,0,2,3))
            t8P = t8 - t8.transpose((1,0,2,3))

            # update t2 amplitudes
            t = t1 + 0.5*t2 + 0.5*t3 + t4P - 0.5*t5P - 0.5*t6P + 0.25*t7 + t8P
            t *= Ep 

            # evaluate Ecorr
            Ec = 0.25 * np.sum(g[o,o,v,v] * t )
            dE = np.fabs(Ec - self.Ec)

            # save stuff
            self.Ec = Ec
            self.t = t

            print       ("@CCD {:3d} {:20.15f} {:20.15f}"  .format(k,Ec,dE) )

            if dE < self.conv: break

        return self.Ec
