import psi4
import numpy as np
import scipy.linalg as la
import itertools as it
import sys
sys.path.append("../../5/avcopan/")
from uhf import UHF


class CIS:

    def __init__(self, mol, mints):
      uhf = UHF(mol, mints)
      uhf.compute_energy()
      self.nocc = uhf.nocc
      self.norb = uhf.norb
      self.ndet = self.nocc * (self.norb - self.nocc)
      self.E0   = uhf.E
      self.e    = uhf.e
      self.g    = transform_tei(uhf.g, uhf.C) # antisymmetrized two-electron integrals, spin-orbital (MO) basis

    def singles_iterator(self):
      return enumerate(it.product(range(0,self.nocc), range(self.nocc,self.norb)))

    def compute_excitation_energies(self):
      ndet, E0, e, g = self.ndet, self.E0, self.e, self.g
      Hn = np.zeros((ndet, ndet))
      for P, (i,a) in self.singles_iterator():
        for Q, (j,b) in self.singles_iterator():
          Hn[P,Q] = (e[a] - e[i])*(i==j)*(a==b) + g[a,j,i,b]
      E, C = la.eigh(Hn)

      print('{:>6s}{:>20s}'.format('State','Energy (Eh)'))
      for K in range(ndet):
        print('{:6d}{:20.15f}'.format(K, E[K]))
      self.E, self.C = E, C
      return E


def transform_tei(gao, C):
  return np.einsum('Pp,Pqrs->pqrs', C, 
           np.einsum('Qq,PQrs->Pqrs', C,
             np.einsum('Rr,PQRs->PQrs', C,
               np.einsum('Ss,PQRS->PQRs', C, gao)
             )
           )
         )
