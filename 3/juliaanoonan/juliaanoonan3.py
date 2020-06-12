#Hartree-Fock, yeah baby

import numpy as np
import psi4.core
import scipy.linalg as la

class rhf(object):

    def __init__(self, mol, mints):

        self.mol = mol
        self.mints = mints

        # initialization 1: read in Vnuc, integrals

        self.Vnuc = mol.nuclear_repulsion_energy()
        self.natom = mol.natom()
        self.charge = mol.molecular_charge() # 0 if neutral, -1 if single anion, etc
        self.norb = mints.basisset().nbf()  # num (spatial) orbitals = num basis functions
        self.nuclear_charges = [mol.Z(A) for A in range(self.natom)] # nuc charges by atom
        self.nelec = sum(self.nuclear_charges) - self.charge #num e- = nuc charges - total charge
        self.nocc = int(self.nelec/2)

        self.S = np.array(mints.ao_overlap()) # overlap integrals
        self.T = np.array(mints.ao_kinetic()) # kinetic energy integral
        self.V = np.array(mints.ao_potential()) # electron-nuclear attraction integrals
        self.g = np.array(mints.ao_eri()).transpose((0,2,1,3)) # electron-electron repulsion integrals, in physicist's notation

        # initialization 2: form orthogonalizer X

        self.X = mints.ao_overlap()
        self.X.power(-0.5, 1.e-16)
        self.X = np.array(self.X)

        # initialization 3: set D = 0 as starting guess

        self.D = np.zeros((self.norb,self.norb))

        # form core Hamiltonian (is not altered)


        self.h = self.T + self.V

        # form nu 

        self.J = np.einsum('abcd, db -> ac', self.g, self.D)         #first part, Coulomb, summing  SR RATHER THAN RS?
        self.K = np.einsum('abdc, db -> ac', self.g, self.D)         #second part, exchange, summing MADE TINY DIFFERENCE
        self.G = self.J - (0.5*self.K)

        # build Fock matrix

        self.f = self.h + self.G

        self.E = 0  #initial energy guess, ridiculous but fine bc it won't converge first try

        self.iter = 0

        self.Convergence = False

    def toconverge(self):  #could easily add parameter of level of convergence

        while self.Convergence != True:

            self.oldf = np.copy(self.f)
            self.oldD = np.copy(self.D)
            self.oldE = np.copy(self.E)

            # remake Fock matrix with new density

            self.J = np.einsum('abcd, db -> ac', self.g, self.D)         #first part, Coulomb, summing
            self.K = np.einsum('abdc, db -> ac', self.g, self.D)         #second part, exchange, summing
            self.G = self.J - (0.5*self.K)

            self.f = self.h + self.G

            # orthogonalize Fock matrix to AO basis

            self.fao = np.matmul(np.matmul(self.X, self.f), self.X)

            # get orbital energies and MO coefficients in terms of eigvals&vecs

            self.eaos, self.caos = np.linalg.eigh(self.fao)

            # backtransform coeffecients to AO basis

            self.C = np.matmul(self.X, self.caos)

            # build density matrix


            self.touse = self.C[:,:self.nocc]  #important! only sum over occupied orbitals
            self.D = np.einsum('pi,qi -> pq', self.touse, np.conj(self.touse)) * 2

            # evaluate energy


            self.summ = self.h + (0.5*self.G)
            self.ee = np.einsum('pq,qp', self.summ, self.D)             #electonic energy
            self.E = self.ee + self.Vnuc
            print(self.E)

            if self.iter > 1:  # so that it has at least two cycles to get away from zero
               if abs(self.E - self.oldE) < 0.0000000001 and np.allclose(self.oldD, self.D, rtol=1e-04, atol=1e-07, equal_nan=True):
                   self.Convergence = True

            self.iter += 1

        print("Energy of {:f} Hartrees, converged in {:d} iterations.".format(self.E, self.iter))
