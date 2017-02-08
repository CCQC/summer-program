#!/usr/bin/env python3

import psi4, configparser, numpy as np
import sys
sys.path.insert(0,"../../3/aewiens")
from rhf import RHF


class MP2:

	def __init__(self,mol,mints):

		rhf       = RHF(mol,mints)
		self.E0   = rhf.computeEnergy()
		self.docc = rhf.docc
		self.nbf  = len(rhf.S) 
		self.e    = rhf.e
		self.G    = rhf.G
		self.C    = rhf.C


	def get_energy(self):

		docc = self.docc
		nbf  = self.nbf
		E0   = self.E0
		e    = self.e

		Gmo   = self.transform_integrals()
		Ecorr = 0.0

		for i in range(docc):
			for j in range(docc):
				for a in range(docc,nbf):
					for b in range(docc,nbf):
						Ecorr += (2*Gmo[i,j,a,b]-Gmo[i,j,b,a])*Gmo[i,j,a,b]/(e[i]+e[j]-e[a]-e[b])

		return E0 + Ecorr


	def transform_integrals(self):
		G,C = self.G, self.C
		return np.einsum("Pqrs,Pp->pqrs",
			np.einsum("PQrs,Qq->Pqrs",
			np.einsum("PQRs,Rr->PQrs",
			np.einsum("PQRS,Ss->PQRs", G, C), C) ,C) ,C)


if __name__ == '__main__':
	config = configparser.ConfigParser()
	config.read('Options.ini')

	molecule   = psi4.geometry( config['DEFAULT']['molecule'] )
	molecule.update_geometry()

	basis = psi4.core.BasisSet.build(molecule, "BASIS", config['DEFAULT']['basis'],puream=0)
	mints = psi4.core.MintsHelper(basis)
	mp2 = MP2(molecule,mints)
	print( mp2.get_energy() )


"""
##  Slower but instructive ways to transform integrals

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
