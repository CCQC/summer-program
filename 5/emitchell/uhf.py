#!/bin/python

import sys
import numpy as np
from numpy import linalg as LA
from scipy import linalg as SLA

np.set_printoptions(threshold=sys.maxsize)

class UHF(object):
        def __init__(self,mol,mints):
                self.mol = mol
                self.mints = mints
                self.Vnu = mol.nuclear_repulsion_energy()
                S = np.array(mints.ao_overlap())    # overlap integrals
                T = np.array(mints.ao_kinetic())    # kinetic energy integrals
                V = np.array(mints.ao_potential())  # electron-nuclear attraction integrals
                g = np.array(mints.ao_eri())   # electron-electron repulsion integrals
                self.S = SLA.block_diag(S,S)   # Converts S to spin-AO basis
                self.T = SLA.block_diag(T,T)   # Converts T to spin-AO basis
                self.V = SLA.block_diag(V,V)   # Converts V to spin-AO basis
                self.G = np.transpose(self.blockTensor(g),(0,2,1,3))   # Converts to spin-AO basis in physicist's notation
                self.X = SLA.fractional_matrix_power(self.S,-0.5) 
                self.nocc = int(sum([mol.Z(A) for A in range(mol.natom())] )) - mol.molecular_charge()   # Number of occupied orbitals
                self.norb = mints.basisset().nbf()   # number of (spatial) orbitals = number of basis functions
                self.h = self.T + self.V   # Core Hamiltonian
                # Guess density matrix
                h = T + V
                self.D = np.kron(h,np.array([[1,0],[0,-1]]))

        def blockTensor(self, G):
                for i in range(2): 
                        G = np.transpose(G,(2,3,0,1))
                        G = np.kron(np.eye(2),G)
                G = np.array(G)
                return G

        def scf(self):
                Dold = np.copy(self.D)
                # Build Fock Matrix
                J = np.einsum('mrns,sr -> mn', self.G, Dold)
                K = np.einsum('mrsn,sr -> mn', self.G, Dold)
                self.v = J - K 
                f = self.h + self.v
                # Diagonalize Fock Matrix
                self.F = LA.multi_dot([self.X, f, self.X])   # Transformed to orthogonalized AO basis
                self.energies, coeffs = LA.eigh(self.F)  # Returns orbital energies and MO coefficients
                # Build Density Matrix
                self.C = np.dot(self.X,coeffs)
                Cocc = self.C[:,:self.nocc]   # Removes virtual orbitals
                self.co = np.copy(Cocc)
                self.D = np.einsum('mi,ni -> mn', Cocc, np.conj(Cocc)) 
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
                print("UHF Final Energy:  {:> 20.14f}".format(Enew))
                return Enew
