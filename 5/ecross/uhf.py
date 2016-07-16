#!/anaconda/bin/python

import psi4 
import numpy as np
import re
from scipy import linalg as la
from numpy import einsum as es

class UHF(object):
    def __init__(self, mol, mints):
        self.mol, self.mints = mol, mints
        self.nbf             = self.mints.basisset().nbf()
        self.maxiter         = psi4.get_global_option('MAXITER')
        self.e_convergence   = psi4.get_global_option('E_CONVERGENCE')
        self.vnu             = mol.nuclear_repulsion_energy()

        # Determine number of electrons
        self.nelec           = -mol.molecular_charge()
        for i in range(mol.natom()):
            self.nelec += int(mol.Z(i))

        # Import one and two electron integrals in spatial AO basis
        self.s_spat = np.array(mints.ao_overlap())
        self.t_spat = np.array(mints.ao_kinetic())
        self.v_spat = np.array(mints.ao_potential())
        self.g_spat = np.array(mints.ao_eri())

        # Transform 1- and 2-e integrals to spin AO basis
        self.s, self.t, self.v, self.g = self.spin_integrals()

        # Form orthogonalizer, X = S^(-0.5)
        self.x = la.sqrtm(la.inv(self.s))

        # Initial HF parameters: D=0, E=0 
        self.d     = np.zeros((2*self.nbf,2*self.nbf))
        self.E     = 0
        self.diff  = 0
        self.count = 0

        # Iterate to self-consistency
        self.iterate()

        # Compare with Psi4 uhf energy
        self.check_answer()

    def spin_integrals(self):
        # Import relevant variables
        s_spat, t_spat = self.s_spat, self.t_spat
        v_spat, g_spat = self.v_spat, self.g_spat
        nbf = self.nbf

        # Block out 2*nbf x 2*nbf (x 2*nbf x 2*nbv) tensors
        s = np.zeros((2*nbf,2*nbf))
        t = np.zeros((2*nbf,2*nbf))
        v = np.zeros((2*nbf,2*nbf))
        g = np.zeros((2*nbf,2*nbf,2*nbf,2*nbf))

        # Transform one electron integrals 
        s[:nbf,:nbf] = s_spat[:nbf,:nbf] 
        s[nbf:,nbf:] = s_spat[:nbf,:nbf]
        t[:nbf,:nbf] = t_spat[:nbf,:nbf] 
        t[:nbf,:nbf] = t_spat[:nbf,:nbf] 
        v[nbf:,nbf:] = v_spat[:nbf,:nbf]
        v[nbf:,nbf:] = v_spat[:nbf,:nbf]

        # Transform two electron integrals
        for i in range(nbf):
            for j in range(nbf):
                # Top left quarter of metamatrix
                g[i,j,nbf:,nbf:]         = g_spat[i,j,:nbf,:nbf]
                g[i,j,:nbf,:nbf]         = g_spat[i,j,:nbf,:nbf]
                # Bottom right quarter of metamatrix
                g[nbf+i,nbf+j,nbf:,nbf:] = g_spat[i,j,:nbf,:nbf]
                g[nbf+i,nbf+j,:nbf,:nbf] = g_spat[i,j,:nbf,:nbf]
        return(s,t,v,g)

    def iterate(self):
        # Import relevant variables
        s, g, d, x        = self.s, self.g, self.d, self.x
        nelec, vnu, count = self.nelec, self.vnu, self.count
        maxiter, e_conv   = self.maxiter, self.e_convergence

        h = self.t + self.v 

        for iterations in range(maxiter):
            count += 1 
            J = es('ijkl,lk->ij',g,d)   #Coulomb integral 
            K = es('iljk,lk->ij',g,d)   #Exchange integral
            f = h + J - K               #Fock matrix
            ft = np.dot(x,np.dot(f,x))  #Transform Fock matrix
            e_orb,ct = la.eigh(ft)      #Diagonalize transformed Fock
            c = np.dot(x,ct)            #Backtransform coef matrix
            co = c[:,:nelec]            #Truncate coef matrix
            cc = np.conj(co)            #Conjugate coef matrix
            d = es('ik,jk->ij',co,cc)   #Build density matrix
            op = h + 0.5*J - 0.5*K      #Construct energy operator
            ee = es('ij,ji',op,d)       #Evaluate electronic energy
            e = vnu + ee                #Evaluate total energy
            diff = abs(self.E - e)
            
            if count == maxiter:
                self.print_failure()
                break

            elif diff < e_conv:
                self.print_success()
                break

            else:
                self.count, self.E, self.d, self.diff = count, e, d, diff
                self.print_iteration()

    def print_failure(self):
        string  = '\n' + '='*78 + '\n'*2
        string += 'Maximum number of iterations reached.\nUnconverged energy: {:16.10f}\tConvergence: {:16.10f}\n'.format(self.E,self.diff)
        string += '\n' + '='*78
        print(string)

    def print_success(self):
        print('\n' + '='*78)
        print('Converged after {} iterations.'.format(self.count))
        print('Final RHF-SCF Energy: {:.10f}'.format(self.E))
        print('='*78)

    def print_iteration(self):
        print('Iteration {:d}\tEnergy: {:16.10f}\tConvergence: {:17.10f}'.format(self.count,self.E,self.diff))

    def check_answer(self):
        output = open('output.dat','r').read()
        psi = float(re.findall('\s+Total Energy\s+\=\s+(\-\d+\.\d+)',output)[0])
        print('\n' + 'Discrepancy with Psi4:\t' + str(psi-self.E))
        
