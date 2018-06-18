#!/anaconda/bin/python

import psi4 
import numpy as np
import re
from scipy import linalg as la
from numpy import einsum as es
from timeit import default_timer

class UHF(object):
    def __init__(self, mol, mints):
        """
        Optimize a spin-orbital basis by Unrestricted Hartree-Fock theory.

        Inputs: 
        1) Molecule object
        2) Integral package

        IMPORTANT MEMBER VARIABLES:
        Optimised energy................self.E
        Optimised density matrix........self.d

        n.b. Chemist's notation is used for in the evaluation of two electron
        integrals.
        """
        self.mol, self.mints = mol, mints
        self.nbf             = self.mints.basisset().nbf()
        self.maxiter         = psi4.get_global_option('MAXITER')
        self.e_convergence   = psi4.get_global_option('E_CONVERGENCE')
        self.vnu             = mol.nuclear_repulsion_energy()

        self.duration = 0
        self.print_header()

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

        n = self.nbf

        # Form orthogonalizer, X = S^(-0.5)
        self.x = la.sqrtm(la.inv(self.s))

        # Initial HF parameters: D=0, E=0 
        self.d     = np.zeros((2*self.nbf,2*self.nbf))
        self.E     = 0
        self.diff  = 0
        self.count = 0

        # Iterate to self-consistency
        self.iterate()

        # The following function compares output to Psi4 energy
        """
        self.check_answer()
        """

    def spin_integrals(self):
        """
        Transform integrals from spatial to spin-orbital basis.

        Outputs: Overlap (s), one electron (t,v), and two electron (g) integrals in
                 spin-orbital AO basis
        """
        s_spat, t_spat = self.s_spat, self.t_spat
        v_spat, g_spat = self.v_spat, self.g_spat
        nbf = self.nbf

        print('Transforming integrals to spin-orbital AO basis...')
        start = default_timer()

        # Block out 2*nbf x 2*nbf (x 2*nbf x 2*nbv) tensors
        s = np.zeros((2*nbf,2*nbf))
        t = np.zeros((2*nbf,2*nbf))
        v = np.zeros((2*nbf,2*nbf))
        g = np.zeros((2*nbf,2*nbf,2*nbf,2*nbf))

        # Transform one electron integrals 
        s[:nbf,:nbf] = s_spat[:nbf,:nbf] 
        s[nbf:,nbf:] = s_spat[:nbf,:nbf]
        t[:nbf,:nbf] = t_spat[:nbf,:nbf] 
        t[nbf:,nbf:] = t_spat[:nbf,:nbf] 
        v[:nbf,:nbf] = v_spat[:nbf,:nbf]
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

        duration = default_timer() - start
        self.duration += duration
        print(' '*40 + '...completed. Runtime: {:s}'.format(str(duration)) + '\n')
        print('='*78)

        return(s,t,v,g)

    def iterate(self):
        """
        Iterate to self consistency.

        After each iteration, the following member variables are updated:
            self.count, self.E, self.d, self.diff
        """
        s, g, d, x        = self.s, self.g, self.d, self.x
        nelec, vnu, count = self.nelec, self.vnu, self.count
        maxiter, e_conv   = self.maxiter, self.e_convergence

        start = default_timer()

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
                duration = default_timer() - start
                self.duration += duration
                print('\n' + 'Runtime: {:s}'.format(str(duration)))
                self.print_success()
                break

            else:
                self.count, self.E, self.d, self.diff = count, e, d, diff
                self.print_iteration()

    def print_header(self):
        print('\n' + '-'*78 + '\n')
        print(' '*27 + 'Unrestricted Hartree-Fock' + '\n')
        print(' '*30 + 'by Elliot Rossomme' + '\n')
        print('-'*78 + '\n')

    def print_output(self):
        print('\n' + '='*78)
        print('{:20s}{:15.10f}'.format('SCF Energy:',self.rhfe))
        print('-'*35)
        print('{:20s}{:15.10f}'.format('FINAL ENERGY:',self.E))
        print('-'*35)
        print('='*78)

    def print_failure(self):
        string  = '\n' + '='*78 + '\n'*2
        string += 'Maximum number of iterations reached.\nUnconverged energy: {:16.10f}\tConvergence: {:16.10f}\n'.format(self.E,self.diff)
        string += '\n' + '='*78
        print(string)

    def print_success(self):
        print('\n' + '='*78)
        print('Converged after {} iterations.'.format(self.count))
        print('Total Runtime: {:s} seconds'.format(str(self.duration)))
        print('Final RHF-SCF Energy: {:.10f}'.format(self.E))
        print('='*78)

    def print_iteration(self):
        print('Iteration {:d}\tEnergy: {:16.10f}\tConvergence: {:17.10f}'.format(self.count,self.E,self.diff))

    def check_answer(self):
        output = open('output.dat','r').read()
        psi = float(re.findall('\s+Total Energy\s+\=\s+(\-\d+\.\d+)',output)[0])
        print('\n' + 'Discrepancy with Psi4:\t' + str(psi-self.E))
        
