#!/bin/python

import sys
import numpy as np
from numpy import linalg as LA
from scipy import linalg as SLA

class CIS(object):
        def __init__(self,uhf):
                self.C = uhf.C   # MO coefficients from Fock diagonalization
                self.energies = uhf.energies   # Energies from Fock diagonalization
                self.G = uhf.G   # Electron-electron repulsion integrals
                self.nbf = 2 * uhf.norb   # Number of (spatial) orbitals = number of basis functions
                self.nocc = uhf.nocc   # Number of occupied orbitals
                self.nvir = self.nbf - self.nocc   # Number of virtual orbitals

        def ao2mo(self):
                C_conj = np.conj(self.C)   # Complex conjugate of the MO coefficients
                ijab = np.einsum('mnrs,mi,nj,ra,sb -> ijab', self.G, C_conj, C_conj, self.C, self.C)   # Creates the integral in the MO basis
                return ijab

        def excite(self):
                self.ex = [(i,a) for i in range(self.nocc) for a in range(self.nocc,self.nbf)]   # Creates a list of excitations possible

        def shift_matrix(self):
                self.excite()
                ajib = np.transpose(self.ao2mo(),(2,1,0,3))   # Rearranges MO integral to correct indices
                H = np.zeros((self.nocc * self.nvir, self.nocc * self.nvir))   # Creates the empty Hamiltonian matrix
                for r, (i,a) in enumerate(self.ex):   # Loops through the excitations to create the elements of the Hamiltonian
                        for s, (j,b) in enumerate(self.ex):
                                H[r,s] = ((self.energies[a] - self.energies[j]) * (i==j) * (a==b)) + (ajib[a,j,i,b] - ajib[a,j,b,i])
                return H

        def ex_contrib(self,mat,limit):
                Ck = np.transpose(mat)   # Converts to coumn vectors from row vectors
                weights = []
                for i, j in enumerate(Ck):
                        string = ""
                        for k, l in enumerate(j):   # Loops through the column vector to get the weights from each excitation
                                val = (l ** 2) * 100
                                if val >= limit:
                                        string += " {:3} ->{:3}  {:5.1f}% ".format(self.ex[k][0],self.ex[k][1],val)
                        weights.append(string)
                return weights

        def get_energy(self,verbose,lim):
                H = self.shift_matrix()
                Ek, coeffs = LA.eigh(H)   # Solves the Hamiltonian giving back the energies of each excitation and the expansion coefficients
                w = self.ex_contrib(coeffs,lim)
                if verbose == True:
                        print("  {:5}  | {:^10}  |  {:s} (contributions above {:2.1f}%)".format("State","Energy","Excitations",lim))
                        print("-" * 72)
                        for i in range(len(Ek)):
                                print("  {:5}  | {: 10.5f}  |{:s}".format(i+1,Ek[i],w[i]))
                        print("-" * 72)
                return Ek
