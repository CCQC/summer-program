#!/anaconda/bin/python

import psi4
import numpy as np
from scipy import linalg as la

class RHF(object):
    def __init__(self, mol, mints):
        """
        Compute the energy of a system by the Restricted Hartree-Fock
        Self-Consistent Field Method.
        =============================================================
      
        An initial guess of D = 0 is employed.
        """
        self.mol = mol
        self.mints = mints
        self.nbf = self.mints.basisset().nbf()
        self.maxiter = psi4.get_global_option('MAXITER')
        self.e_convergence = psi4.get_global_option('E_CONVERGENCE')

        #Import necessary integrals from psi4
        self.Vnu = mol.nuclear_repulsion_energy()
        self.S = np.matrix(mints.ao_overlap())
        self.T = np.matrix(mints.ao_kinetic())
        self.V = np.matrix(mints.ao_potential())
        self.g = np.array(mints.ao_eri()).transpose((0,2,1,3)) #Import eri integrals in physicist's not'n

        #Construct orthogonalizer, X = S^(-0.5)
        self.X = la.inv(la.sqrtm(self.S))

        #Initial HF guess with D = 0
        self.D = np.zeros((self.nbf,self.nbf))

        # Determine the number of electrons and the number of doubly occupied orbitals
        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        if mol.multiplicity() != 1 or self.nelec % 2:
            raise Exception("This code only allows closed-shell molecules.")
        self.ndocc = self.nelec / 2

        self.E = 0
        self.C = 0
        self.orb_e = 0
        self.count = 0 
        self.diff = 0
        self.iterate()

    def iterate(self):
        """
        Input: self

        Operations: Build Fock matrix (f), transform f to orthogonal AO basis (ft), diagonalize ft,
        transform eigenvectors of ft (ct) to original AO basis (c), truncate c (co) build density
        matrix (d), evaluate energy (e), evaluate convergence (diff).

        Repeat above process until diff is less than specified convergence.

        Output: self.E = converged energy, self.D = final density matrix 
        """
        self.print_header()

        diff = 1
        count = 0
        while diff > self.e_convergence:
            for i in range(self.maxiter):
                count +=1 
                J = np.einsum('upvs,ps->uv',self.g,self.D)    #Calculate Coulomb integrals    output: np array
                K = np.einsum('upsv,ps->uv',self.g,self.D)    #Calculate Exchange integrals   output: np array
                f = self.T + self.V + J - 0.5*K               #Construct Fock matrix          output: np matrix
                ft = self.X * f * self.X                      #Transform Fock matrix          output: np matrix
                nrg,ct_mat = la.eigh(ft)                      #Diagonalize transformed Fock   output: np array, np array
                c = np.dot(self.X,ct_mat)                     #Backtransform coef matrix      output: np array
                co = np.matrix(c[:,:self.ndocc])              #Truncate coef matrix           output: np matrix
                d = 2*np.einsum('ik,jk->ij',co,np.conj(co))   #Construct density matrix       output: np array
                H = 0.5 * (self.T + self.V + f)               #Construct HF Hamiltonian       output: np matrix
                e = self.Vnu + np.einsum('ij,ji',H,d)         #Calculate e from Vnu, H, D     output: float
                diff = abs(self.E - e)                        #Calculate convergence          output: float 
                
                if count == self.maxiter:
                    #string = '\n' + '='*78 + '\n'*2
                    #string += 'Maximum number of iterations reached.\nUnconverged Energy: {:16.10f}\tConvergence: {:16.10f}\n'.format(e,diff)
                    #string += '\n' + '='*78
                    #raise Exception(string)
                    self.print_failure()
                    break
                elif diff < self.e_convergence:
                    self.print_success()
                    break
                else:
                    self.D, self.E, self.C, self.orb_e, self.count, self.diff = d, e, c, nrg, count, diff
                    self.print_iteration()

    def print_header(self):
        print('\n' + '-'*78 + '\n')
        print(' '*10 + 'Restricted Hartree-Fock Self-Consistent Field Optimization' + '\n')
        print(' '*30 + 'by Elliot Rossomme' + '\n')
        print('-'*78 + '\n')

    def print_iteration(self):
        print('Iteration {:d}\tEnergy: {:16.10f}\tConvergence: {:17.10f}'.format(self.count,self.E,self.diff))

    def print_failure(self):
        string = '\n' + '='*78 + '\n'*2
        string += 'Maximum number of iterations reached.\nUnconverged Energy: {:16.10f}\tConvergence: {:16.10f}\n'.format(self.E,self.diff)
        string += '\n' + '='*78
        print(string)

    def print_success(self):
        print('\n' + '='*78)
        print('Converged after {} iterations.'.format(self.count))
        print('Final RHF-SCF Energy: {:.10f}'.format(self.E))
        print('='*78)

        
