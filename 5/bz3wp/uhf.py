# Unrestricted Hartree-Fock class using block diagonalized spin orbital functions

import psi4
import numpy as np
import configparser 
import scipy.linalg as la

class UHF:
    def __init__(self, filename = 'Options.ini'):
        
        config = configparser.ConfigParser()                    # configparser to read in input options
        config.read(filename)
        self.mol = psi4.geometry(config['DEFAULT']['molecule'])
        self.mol.update_geometry()
        basis = psi4.core.BasisSet.build(self.mol,'BASIS', config['DEFAULT']['basis'])
        mints = psi4.core.MintsHelper(basis)                    # mints from psi4 for calculating integrals
        
        self.max_iter = int(config['SCF']['max_iter'])
        self.nelec = -self.mol.molecular_charge()
        for A in range(self.mol.natom()):
            self.nelec += int(self.mol.Z(A))
        self.nocc = self.nelec                                  # number of occupied orbitals
        self.nto = 2* mints.basisset().nbf()                    # total number of spin orbitals
        
        # using blocking functions defined below
        S = block_oei(mints.ao_overlap())                       # blocked one-electron overlap integral
        T = block_oei(mints.ao_kinetic().to_array())            # blocked one-electron kinetic energy integral
        V = block_oei(mints.ao_potential().to_array())          # blocked one-electron potential energy integral
        G = block_tei(mints.ao_eri().to_array())                # blocked two-electron integral
        
        
        self.H = T + V                                          # hamiltonian
        self.g = G.transpose(0,2,1,3) - G.transpose(0,2,3,1)    # antisymmetrized two-electron integrals in physicists' notation
        self.A = np.matrix(la.inv(la.sqrtm(S)))                 # orthogonalizer

        self.C = np.zeros_like(self.H)                          # MO coefficients
        self.e = np.zeros(len(self.H))                          # orbital energies
        self.E_SCF = 0.0
        

    def get_energy(self):
        mol, max_iter, g, H, A, C, e, nocc = self.mol, self.max_iter,self.g, self.H, self.A, self.C, self.e, self.nocc
        
        E_old = 0.0                                             # defining initial values
        D_old = np.zeros_like(H)
        
        for iteration in range(1, self.max_iter+1):             # starting SCF iterations
            
            Gao = np.einsum('pqrs,sq', g, D_old)                # combined coulomb and exchange
            F = H + Gao                                         # constructing fock matrix

            Ft = A.dot(F).dot(A)                                # transform fock matrix
            e, C = np.linalg.eigh(Ft)                           # diagonalize
            C = A.dot(C)                                        # form new SCF eigenvector
            Cnocc = C[:,:nocc]
            D = np.einsum('pi, qi->pq', Cnocc, Cnocc)           # form new density matrix
            
            E_SCF = (1/2)*(np.einsum('pq, pq->', F+H, D)) + mol.nuclear_repulsion_energy()  #calculate SCF energy
            print('UHF iteration {:3d}: energy {:20.14f} dE {:1.5E}'.format(iteration, E_SCF, (E_SCF - E_old)))

            if (abs(E_SCF - E_old) < 1.e-10):                   # convergence criteria
                break
             
            E_old = E_SCF
            D_old = D

        self.e = e                                              # saving values
        self.C = C
        self.E_SCF = E_SCF


# spin blocking functions to transform spatial orbitals to spin orbitals
def block_oei(A):                   # block one-electron integrals
    A = la.block_diag(A, A)
    return A

def block_tei(A):                   # block two-electron integrals
    I = np.identity(2)
    A = np.kron(I, A)
    return np.kron(I, A.T)

# testing 
if __name__=='__main__':
    
    uhf = UHF('Options.ini')
    uhf.get_energy()
