import psi4
import numpy as np

class RMP2(object):

  def __init__(self, rhf):
    G = rhf.g    # AO basis two-electron-integrals
    C = rhf.C    # MO coefficients
    self.e    = rhf.e    # orbital energies
    self.nocc = rhf.nocc # number of occupied orbitals
    self.norb = rhf.norb # total number of orbitals

    # transform integrals to MO basis
    mul = lambda T, M : np.tensordot(T, M, axes=(0,0))
    self.g = mul(mul(mul(mul(G, C), C), C), C)

  def compute_energy(self):
    norb, nocc, e, g = self.norb, self.nocc, self.e, self.g
    E = 0.0
    for i in range(nocc):
      for j in range(nocc):
        for a in range(nocc, norb):
          for b in range(nocc, norb):
            E += g[i,j,a,b] * ( 2*g[i,j,a,b] - g[i,j,b,a] ) / (e[i] + e[j] - e[a] - e[b])

    psi4.print_out('@RMP2 correlation energy: {:20.15f}\n'.format(E))
    self.E = E
    return E

