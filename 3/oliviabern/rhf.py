#!/usr/bin/python

import psi4
import numpy as np
import scipy.linalg as la


class RHF:
    """
    Restricted Hartree-Fock class for obtaining the restricted Hartree-Fock
    energy
    """

    def __init__(self, mol, mints):
        """
        Initialize the rhf
        :param mol: a Psi4 molecule object
        :param mints: a molecular integrals object (from MintsHelper)
        """
        self.mol = mol
        self.mints = mints

        self.V_nuc = mol.nuclear_repulsion_energy()
        self.T = np.matrix(mints.ao_kinetic())
        self.S = np.matrix(mints.ao_overlap())
        self.V = np.matrix(mints.ao_potential())
    

        self.g = np.array(mints.ao_eri()).transpose((0,2,1,3))

        # Determine the number of electrons and the number of doubly occupied orbitals
        self.nelec = -mol.molecular_charge() #accounts for # of e- in ions
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        if mol.multiplicity() != 1 or self.nelec % 2:
            raise Exception("This code only allows closed-shell molecules")
        self.ndocc = self.nelec / 2 #number of doubly occupied orbitals

        self.maxiter = psi4.get_global_option('MAXITER')
        self.e_convergence = psi4.get_global_option('E_CONVERGENCE')

        self.nbf = mints.basisset().nbf()
        self.D = np.zeros((self.nbf,self.nbf))

        #Form orthogonalizer X
        e_values, e_vectors = la.eig(self.S)
        v = np.diagflat(e_values)
        vsqrt = np.sqrt(v)
        vect_inv = la.inv(e_vectors)
        b = np.dot(e_vectors, vsqrt)
        Ssqrt = np.dot(b,vect_inv)
        X = la.inv(Ssqrt)
        self.X = X.real


    def compute_energy(self):
        """
        Compute the rhf energy
        :return: energy
        """
        for i in range(self.maxiter):
            h = self.T + self.V
            j = np.einsum('mrns,rs',self.g,self.D)
            k = np.einsum('msrn,rs',self.g,self.D)
            v = j-.5*k
            f = h + v
            ft = np.dot(self.X , np.dot(f, self.X))
            e , Ct = la.eigh(ft)
            C = np.dot(self.X,Ct)
            OC = C[:,:self.ndocc]
            self.D = 2*np.dot(OC, OC.T)
            T = h + .5*v
            E = np.dot(T,self.D)
            print(np.trace(E)+self.V_nuc)
            
           




