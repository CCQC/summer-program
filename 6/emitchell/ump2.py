#!/bin/python

import sys
import numpy as np
from numpy import linalg as LA
from scipy import linalg as SLA

class UMP2(object):
        def __init__(self,uhf):
                self.C = uhf.C   # MO coefficients from Fock diagonalization
                self.energies = uhf.energies   # Energies from Fock diagonalization
                self.G = uhf.G   # Electron-electron repulsion integrals
                self.nbf = 2 * uhf.norb   # Number of (spatial) orbitals = number of basis functions
                self.nocc = uhf.nocc   # Number of occupied orbitals
                self.nvir = self.nbf - self.nocc   # Number of virtual orbitals

        def ao2mo(self):
                Cocc_conj = np.conj(self.C[:,:self.nocc])   # Complex conjugate of the occupied MO coefficients
                Cvir = self.C[:,self.nocc:]   # Virtual MO coefficients
                ijab = np.einsum('mnrs,mi,nj,ra,sb -> ijab', self.G, Cocc_conj, Cocc_conj, Cvir, Cvir)   # Creates the integral representing double excitations
                return ijab

        def ump2_energy(self):
                ijab = self.ao2mo()
                abij = np.transpose(ijab,(2,3,0,1))   # Conjugate of the MO integral
                Eocc = self.energies[:self.nocc]   # Orbital energies of occupied orbitals
                Evir = self.energies[self.nocc:]   # Orbital energies of virtual orbitals
                E_2 = 0
                for i in range(self.nocc):   # Nested loop represents all possible double excitations
                        for j in range(self.nocc): 
                                for a in range(self.nvir):
                                        for b in range(self.nvir):
                                                E_2 += 0.5 * (ijab[i,j,a,b] * (abij[a,b,i,j] - abij[b,a,i,j])) / (Eocc[i] + Eocc[j] - Evir[a] - Evir[b])
                return E_2   # Total MP2 Correlation Energy
