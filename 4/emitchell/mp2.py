#!/bin/python

import sys
import numpy as np

class MP2(object):
        def __init__(self,rhf):
                self.C = rhf.C   # MO coefficients from Fock diagonalization
                self.energies = rhf.energies   # Energies from Fock diagonalization
                self.G = rhf.g   # Electron-electron repulsion integrals
                self.nbf = rhf.norb   # Number of (spatial) orbitals = number of basis functions
                self.nocc = rhf.nocc   # Number of occupied orbitals
                self.nvir = self.nbf - self.nocc   # Number of virtual orbitals
        def ao2mo(self):
                Cocc_conj = np.conj(self.C[:,:self.nocc])   # Complex conjugate of the occupied MO coefficients
                Cvir = self.C[:,self.nocc:]   # Virtual MO coefficients
                self.IJAB = np.einsum('mnrs,mI,nJ,rA,sB -> IJAB', self.G, Cocc_conj, Cocc_conj, Cvir, Cvir)   # Creates the integral representing double excitations
        def mp2_energy(self):
                self.ao2mo()
                ABIJ = np.transpose(self.IJAB,(2,3,0,1))   # Conjugate of the MO integral
                Eocc = self.energies[:self.nocc]   # Orbital energies of occupied orbitals
                Evir = self.energies[self.nocc:]   # Orbital energies of virtual orbitals
                E_2 = 0
                for i in range(self.nocc):   # Nested loop represents all possible double excitations
                        for j in range(self.nocc): 
                                for a in range(self.nvir):
                                        for b in range(self.nvir):
                                                E_2 += (self.IJAB[i,j,a,b] * ((2 * ABIJ[a,b,i,j]) - ABIJ[b,a,i,j])) / (Eocc[i] + Eocc[j] - Evir[a] - Evir[b])
                return E_2   # Total MP2 Correlation Energy
