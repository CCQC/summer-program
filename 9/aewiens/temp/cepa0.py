import numpy as np
import sys
import psi4
sys.path.insert(0,"../../5/avcopan/")

class CEPA0:

    def __init__(self,uhf):

        self.conv = psi4.get_global_option('E_CONVERGENCE')
        self.maxiter = uhf.maxiter
        
        G = uhf.g
        C = uhf.C
        self.e = uhf.e
        self.norb = uhf.norb
        self.nocc = uhf.nocc
        self.nvirt = self.norb - self.nocc

    # transform integrals to MO basis
        mul = lambda T, M : np.tensordot(T, M, axes=(0,0))
        self.g = mul(mul(mul(mul(G, C), C), C), C)

        #print( np.linalg.norm(self.g) )

        self.Ec = 0.0 


    def compute_energy(self):
        """
        :Return: CCD correlation energy
        """
        g = self.g
        #g = self.transform_integrals(self.G,self.C)
        nocc, nvirt, e = self.nocc, self.nvirt, self.e
        t = np.zeros((self.nocc,self.nocc,self.nvirt,self.nvirt))

        o = slice(None,nocc)
        v = slice(nocc,None)
        x = np.newaxis
        Ep = 1.0/ (e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v] )
        print( np.linalg.norm(Ep) )

        for i in range(2):
            ##  terms

            ##  debug
            #print( np.linalg.norm(g) )
            print(np.linalg.norm(g[o,o,v,v]) )

            t1 = g[o,o,v,v]
            t2 = np.einsum("abcd,ijcd->ijab",g[v,v,v,v],t)
            t3 = np.einsum("klij,klab->ijab",g[o,o,o,o],t)
            t4 = np.einsum("akic,jkbc->ijab",g[v,o,o,v],t)

            ##  permutations
            t4P = t4 - t4.transpose((1,0,2,3)) - t4.transpose((0,1,3,2)) + t4.transpose((1,0,3,2))

            ##  update t
            t = t1 + 0.5*t2 + 0.5*t3 + t4P

            t*= Ep 
            #print( np.linalg.norm( t) )
          
            # evaluate Ecorr
            Ec = 0.25 * np.sum(g[o,o,v,v] * t )
            dE = np.fabs(Ec - self.Ec)

            # save stuff
            self.Ec = Ec
            self.t = t

            print       ("@CEPA {:3d} {:20.15f} {:20.15f}"  .format(i,Ec,dE) )

            if dE < self.conv: break

        return self.Ec

