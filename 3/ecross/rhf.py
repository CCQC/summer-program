#!/anaconda/bin/python

import psi4
import numpy as np
from scipy import linalg as la

class RHF(object):
    def __init__(self, mol, mints):
        """
        Computes the energy of a system by the Restricted Hartree-Fock
        Self-Consistent Field Method.
        ==============================================================
      
        An initial guess of D = 0 is employed.
        """
        self.mol = mol
        self.mints = mints

        self.size = self.mints.basisset().nbf()
        self.Vnu = mol.nuclear_repulsion_energy()
        self.S = np.matrix(mints.ao_overlap())
        self.T = np.matrix(mints.ao_kinetic())
        self.V = np.matrix(mints.ao_potential())
        self.g = np.array(mints.ao_eri())
       
        self.X = self.orthogonalizer()
        self.D = self.density_initial()
        self.col_ex,self.f = self.build_fock()
        self.ft = self.transform_fock()

        # Determine the number of electrons and the number of doubly occupied orbitals
        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        if mol.multiplicity() != 1 or self.nelec % 2:
            raise Exception("This code only allows closed-shell molecules")
        self.ndocc = self.nelec / 2

        self.nrg,self.ct_mat = self.diagonalize()
        self.c_mat = self.back_transform()
        self.D = self.build_density()
        self.E = self.evaluate_energy()

    def orthogonalizer(self):
        """Construct orthogonalizer (X) from overlap matrix."""
        X = la.inv(la.sqrtm(self.S))
        return(X)

    def density_initial(self):
        """Construct density matrix with initial guess D = 0."""
        D = np.zeros((self.size,self.size))
        return(D)
    
    def build_fock(self):
        col_ex = np.zeros((self.size,self.size))
        for u in range(self.size):
            for v in range(self.size):
                value = 0
                for p in range(self.size):
                    for s in range(self.size):
                        value += self.g[u,v,p,s] - 0.5 * self.g[u,s,p,v]
                        value *= self.D[s,p]
                col_ex[u,v] += value
        f = self.T + self.V + col_ex
        return(col_ex,f)

    def transform_fock(self):
        ft = la.inv(self.X) * self.f * self.X
        return(ft)

    def diagonalize(self):
        nrg,ct_mat = la.eigh(self.ft)
        return(nrg,ct_mat)

    def back_transform(self):
        c_mat = self.X * self.ct_mat 
        return(c_mat) 

    def build_density(self):
        for i in range(self.size):
            for j in range(self.size):
                value = 0
                for k in range(self.ndocc):
                    value += 2 * self.c_mat[i,k] * np.conj(self.c_mat[j,k])
                self.D[i,j] = value
        return(self.D)

    def evaluate_energy(self):
        E = self.Vnu
        for i in range(self.size):
            for j in range(self.size):
                E += (self.T[i,j] + self.V[i,j] + 0.5 * self.col_ex[i,j]) * self.D[j,i]
        return(E)

    def iterate(self):
        e = self.evaluate_energy()
        diff = np.absolute(e - self.E)
        
        for i in range(maxiter):
            if diff < self.e_convergence:
                pass
            else:
                self.build_fock()
                self.transform_fock()
                self.diagonalize()
                self.back_transform()
                self.build_density()
                e = self.evaluate_energy()


        
       

