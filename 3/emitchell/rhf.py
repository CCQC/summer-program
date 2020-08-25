#!/bin/python

import sys
import numpy as np
from numpy import linalg as LA
from scipy import linalg as SLA

class RHF(object):
        def __init__(self,mol,mints):
                self.mol = mol
                self.mints = mints
                self.Vnu = mol.nuclear_repulsion_energy()
                self.S = np.array(mints.ao_overlap())    # overlap integrals
                self.T = np.array(mints.ao_kinetic())    # kinetic energy integrals
                self.V = np.array(mints.ao_potential())  # electron-nuclear attraction integrals
                self.g = np.transpose(np.array(mints.ao_eri()),(0,2,1,3))        # electron-electron repulsion integrals, changin indices to match equations (physicist's notation)
                self.X = SLA.fractional_matrix_power(self.S,-0.5)   # Orthoganlizer
                self.nocc = int(sum([mol.Z(A) for A in range(mol.natom())] ) / 2)   # Number of occupied orbitals
                self.norb = mints.basisset().nbf() # number of (spatial) orbitals = number of basis functions
                self.h = self.T + self.V   # Core Hamiltonian
                # Guess density matrix 
                self.D = np.zeros((self.norb,self.norb))

        def scf(self):
                Dold = np.copy(self.D)
                # Build Fock Matrix
                J = np.einsum('mrns,rs -> mn', self.g, Dold)
                K = np.einsum('mrsn,rs -> mn', self.g, Dold)
                self.v = J - (0.5 * K)
                f = self.h + self.v
                # Diagonalize Fock Matrix
                self.F = LA.multi_dot([self.X, f, self.X])   # Transformed to orthogonalized AO basis
                self.energies, coeffs = LA.eigh(self.F)  # Returns orbital energies and MO coefficients
                self.co = np.copy(coeffs)
                # Build Density Matrix
                self.C = np.dot(self.X,coeffs)
                Cocc = self.C[:,:self.nocc]   # Removes virtual orbitals
                self.D = 2 * np.einsum('mi,ni -> mn', Cocc, np.conj(Cocc)) 
                self.get_energy()
                return self.Etot

        def get_energy(self):   # Determines the electronic energy from density matrix to get total energy
                Ee1 = np.einsum('mn,nm', self.h, self.D)
                Ee2 = np.einsum('mn,nm', self.v, self.D)
                Ee = Ee1 + (0.5 * Ee2)
                self.Etot = Ee + self.Vnu
                return self

        def optimize_energy(self,conv_tol,maxiter,verbose):   # Optimizes the energy with respect to the density matrix
                Eold = 0
                i = 0
                conv = 100
                if verbose == True:
                        print("-" * 49)
                while abs(conv) > float(conv_tol):
                        Enew = self.scf()
                        if verbose == True:
                                print("Iteration: {:4}  |  Energy: {: 20.14f}".format(i,Enew))
                        conv = Eold - Enew
                        Eold = Enew
                        if i == int(maxiter):
                                sys.exit("Could not converge in {} iterations.".format(maxiter)) 
                        i += 1
                if verbose == True:
                        print(("-" * 49))
                        print("RHF Final Energy:  {:> 20.14f}".format(Enew))
                return Enew
