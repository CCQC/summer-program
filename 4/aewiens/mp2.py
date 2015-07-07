import numpy as np
from scipy import linalg as la
from psi4_helper import get_docc, get_nbf, get_conv, get_maxiter

class MP2:
    def __init__(self,mol,mints):
        self.docc = get_docc(mol)
        self.nbf = get_nbf(mints)
        self.conv = get_conv()
        self.maxiter = get_maxiter()
        self.S = np.matrix(mints.ao_overlap() )
        self.X = np.matrix(la.funm(self.S,lambda x : x**(-0.5)))
        self.T = np.matrix(mints.ao_kinetic() )
        self.V = np.matrix(mints.ao_potential() )
        self.H = self.T+self.V
        self.G = np.array(mints.ao_eri() ).swapaxes(1,2)
        self.Vnu = mol.nuclear_repulsion_energy()
        self.D = np.matrix(np.zeros(self.S.shape))
        self.E = 0

    def get_energy(self):

        #rename class variables
        docc, nbf, conv, maxiter = self.docc, self.nbf, self.conv, self.maxiter
        S, X, T, V, H, G, Vnu, D = self.S, self.X, self.T, self.V, self.H, self.G, self.Vnu, self.D
        E = self.E

        # get SCF energy
        for i in range(maxiter):
            J= np.einsum("ijkl,jl->ik",G,D)
            K = np.einsum("ijkl,jk->il",G,D)
            F = H+J-0.5*K      # start with F = Hcore
            tF = X*F*X
            e, tC = la.eigh(tF)
            C = X*tC
            oC = C[:,:docc]
            D = 2*oC*oC.T
            E0 = E
            E = np.trace(0.5*(H+F)*D)+Vnu
            dE = np.fabs(E-E0)
            self.E = E
            #print E
            #print('{:20.15f}{:20.15f}'.format(E,dE))
            if dE < self.conv : break

        Gmo = np.zeros(G.shape)
        for p in range(nbf):
          for q in range(nbf):
            for r in range(nbf):
              for s in range(nbf):
                for i in range(nbf):
                 for j in range(nbf):
                   for k in range(nbf):
                     for l in range(nbf):
                       Gmo[p,q,r,s] += C[i,p]*C[j,q]*C[k,r]*C[l,s]*G[i,j,k,l]

        Ecorr = 0.0
        for i in range(docc):
          for j in range(docc):
            for a in range(docc,nbf):
              for b in range(docc,nbf):
                Ecorr += (2*Gmo[i,j,a,b]-Gmo[i,j,b,a])*Gmo[i,j,a,b]/(e[i]+e[j]-e[a]-e[b])

        return E + Ecorr
