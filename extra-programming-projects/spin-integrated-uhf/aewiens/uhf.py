#!/usr/bin/env python

import numpy as np
from scipy import linalg as la

from psi4_helper import get_maxiter, get_nbf, get_nelec, get_conv


class UHF:

    def __init__(self,mol,mints):

        mult = mol.multiplicity()
        nelec = get_nelec(mol)       
        self.conv = get_conv()

        self.Na = int( 0.5*(nelec+mult-1) )
        self.Nb = nelec - self.Na

        V = np.array( mints.ao_potential() )
        T = np.array( mints.ao_kinetic() )
        G = np.array( mints.ao_eri() )
        self.S = np.array( mints.ao_overlap() )

        self.Hcore = T + V
        self.G = G.transpose((0,2,1,3))
        self.X = np.matrix( la.funm(self.S, lambda x : x**(-0.5) ) )
        self.Vnu = mol.nuclear_repulsion_energy()

        ##  alpha and beta density matrices
        """
        self.Da = np.zeros((self.X.shape))
        self.Db = np.zeros((self.X.shape))
        """
        self.Da = np.random.rand(*self.X.shape)
        self.Da = self.Da + self.Da.T

        self.Db = np.random.rand(*self.X.shape)
        self.Db = self.Db + self.Db.T

        self.E = 0.0

    def compute_energy(self):

        Hcore, G, Da, Db, X = self.Hcore, self.G, self.Da, self.Db, self.X

        for i in range( get_maxiter() ):
            
            va = np.einsum("mnrs,ns->mr",G,Da) - np.einsum("mnsr,ns->mr",G,Da) + np.einsum("mnrs,ns->mr",G,Db)
            vb = np.einsum("mnrs,ns->mr",G,Db) - np.einsum("mnsr,ns->mr",G,Db) + np.einsum("mnrs,ns->mr",G,Da)

            Fa = Hcore + va
            Fb = Hcore + vb

            tFa = X.dot(Fa.dot(X))
            tFb = X.dot(Fb.dot(X))

            ea, tCa = la.eigh(tFa)
            eb, tCb = la.eigh(tFb)

            Ca = X.dot(tCa)
            Cb = X.dot(tCb)

            oCa = Ca[:,:self.Na]
            oCb = Cb[:,:self.Nb]

            Da = oCa.dot(oCa.T)
            Db = oCb.dot(oCb.T)

            E0 = self.E
            E = np.trace( (Hcore+0.5*va).dot(Da) ) + np.trace( (Hcore+0.5*vb).dot(Db) ) + self.Vnu
            dE = np.fabs(E-E0)

            if dE < self.conv: break

            print("UHF  {:>4} {: >21.13}  {: >21.13}".format(i,E,dE))

            ## save
            self.Da, self.Db, self.E = Da, Db, E

        return self.E

