import psi4
import numpy as np
import scipy.linalg as la
import sys
sys.path.append("../../5/avcopan/")
from uhf import UHF


class MP2:

    def __init__(self, mol, mints):
      uhf = UHF(mol, mints)
      uhf.compute_energy()
      self.g   = transform_tei(uhf.g, uhf.C) # antisymmetrized two-electron integrals, spin-orbital (MO) basis
      self.uhf = uhf

    def compute_energy(self):
      uhf = self.uhf
      e, nocc, nbf, g = uhf.e, uhf.nocc, uhf.nbf, self.g
      Ec = 0.0
      for i in range(0,nocc):
        for j in range(0,nocc):
          for a in range(nocc,nbf):
            for b in range(nocc,nbf):
              Ec += (1./4) * g[i,j,a,b]**2 / (e[i]+e[j]-e[a]-e[b])
      self.E = uhf.E + Ec
      print('{:7s}{:>20s}{:>20s}'.format('MP2','E','Ec'))
      print('{:27.15f}{:20.15f}'.format(self.E, Ec))
      return Ec

def transform_tei(gao, C):
  # g_pqrs = sum_P C_Pp (sum_Q C_Qq (sum_R C_Rr (sum_S C_Ss gao_PQRS)))
  return np.einsum('Pp,Pqrs->pqrs', C, 
           np.einsum('Qq,PQrs->Pqrs', C,
             np.einsum('Rr,PQRs->PQrs', C,
               np.einsum('Ss,PQRS->PQRs', C, gao)
             )
           )
         )



# Other ways of doing the integral transformation:

def transform_tei_einsum(gao, C):
  # g_pqrs = sum_PPQRS gao_PQRS C_Pp C_Qq C_Rr C_Ss
  return np.einsum('PQRS,Pp,Qq,Rr,Ss->pqrs', gao, C, C, C, C)

def transform_tei_forloop_n8(gao, C):
  nbf = gao.shape[0]
  g = np.zeros(gao.shape)
  for p in range(nbf):
    print(p)
    for q in range(nbf):
      print(q)
      for r in range(nbf):
        print(r)
        for s in range(nbf):
          for mu in range(nbf):
            for nu in range(nbf):
              for rh in range(nbf):
                for si in range(nbf):
                  g[p,q,r,s] += gao[mu,nu,rh,si] * C[mu,p] * C[nu,q] * C[rh,r] * C[si,s]
  return g

def transform_tei_forloop_n5(gao, C):
  nbf = gao.shape[0]
  g1 = np.zeros(gao.shape)
  g2 = np.zeros(gao.shape)
  for mu in range(nbf):
    for nu in range(nbf):
      for rh in range(nbf):
        for s in range(nbf):
          for si in range(nbf):
            g1[mu,nu,rh,s] += gao[mu,nu,rh,si] * C[si,s]
  for mu in range(nbf):
    for nu in range(nbf):
      for r in range(nbf):
        for s in range(nbf):
          for rh in range(nbf):
            g2[mu,nu,r,s] += g1[mu,nu,rh,s] * C[rh,r]
  g1.fill(0)
  for mu in range(nbf):
    for q in range(nbf):
      for r in range(nbf):
        for s in range(nbf):
          for nu in range(nbf):
            g1[mu,q,r,s] += g2[mu,nu,r,s] * C[nu,q]
  g2.fill(0)
  for p in range(nbf):
    for q in range(nbf):
      for r in range(nbf):
        for s in range(nbf):
          for mu in range(nbf):
            g2[p,q,r,s] += g1[mu,q,r,s] * C[mu,p]
  return g2
