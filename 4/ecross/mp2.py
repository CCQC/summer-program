#!/anaconda/bin/python

import psi4
import numpy as np
from scipy import linalg as la
from timeit import default_timer

# Import RHF class
import sys
sys.path.append('/Users/ecrossomme/git/summer-program/3/ecross/')
from rhf import RHF 


class MP2(object):
    def __init__(self, mol, mints):
        self.mol = mol
        self.mints = mints 

        # Set global variables for occupied and unoccupied orbitals
        self.nbf = self.mints.basisset().nbf()
        
        # Import RHF energy, two electron integrals, coefficient matrix, and orbital energies
        self.rhfe, self.g, self.c, self.orb_e = self.rhf_energy()

        self.nocc, self.nvir, self.occ_e, self.vir_e = self.occupancy()

        self.print_header()

        self.mp2e = self.mp2_energy()
        self.E = self.rhfe + self.mp2e
        self.print_output()
        
    def rhf_energy(self):
        """
        Optimize RHF orbitals; return RHF energy, two-electron integrals,
        and coefficient matrix
        """
        rhf = RHF(self.mol, self.mints)
        return(rhf.E,rhf.g,rhf.C,rhf.orb_e)

    def occupancy(self):
        """
        Determine orbital occupancies; raise exception for open-shell systems.
        """
        nelec = -self.mol.molecular_charge()
        for A in range(self.mol.natom()):
            nelec += int(self.mol.Z(A))
        if self.mol.multiplicity() != 1 or nelec % 2:
            raise Exception('This code only allows closed-shell molecules')
        nocc = nelec / 2
        nvir = self.nbf - nocc
        occ_e = self.orb_e[:nocc]
        vir_e = self.orb_e[nocc:]
        return(nocc,nvir,occ_e,vir_e)

    def mp2_energy(self):
        """
        1) Transform 2 electron integrals from AO to MO basis
        2) Compute MP2 correction to energy
        ---
        Return: MP2 correction to energy
        """
        nocc, nvir, nbf, g = self.nocc, self.nvir, self.nbf, self.g
        c, occ_e, vir_e = self.c, self.occ_e, self.vir_e
        ct = c.transpose()
        c_occ = c[:,:nocc]
        c_vir = c[:,nocc:]

        print('Transforming integrals to MO basis...')
        start = default_timer()

        # Transform 2 electron integrals to MO basis
        # FOUR EINSUM FUNCTIONS; 0.0339229106903 seconds for H2O cc-pVDZ
        gt = np.einsum("ijkl,iu->ujkl",
             np.einsum("ijkl,jv->ivkl",
             np.einsum("ijkl,kp->ijpl",
             np.einsum("ijkl,ls->ijks",g,c_vir),c_vir),c_occ),c_occ)
        duration = str(default_timer() - start)

        # Evaluate mp2 energy
        print(' '*40 + '...completed. Runtime: {:s}'.format(duration))
        print('Evaluating MP2 contribution to energy...')
        start = default_timer()
        mp2e = 0
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nvir):
                    for b in range(nvir):
                        num = gt[i,j,a,b] * (2*gt[i,j,a,b] - gt[i,j,b,a])
                        den = occ_e[i] + occ_e[j] - vir_e[a] - vir_e[b]
                        mp2e += num / den
        duration = str(default_timer() - start)
        print(' '*39 + ' ...completed. Runtime: {:s}'.format(duration))
        return(mp2e)

    
    def print_header(self):
        print('\n' + '-'*78 + '\n')
        print(' '*32 + 'MP2 Correction' + '\n')
        print(' '*30 + 'by Elliot Rossomme' + '\n')
        print('-'*78 + '\n')

    def print_output(self):
        print('\n' + '='*78)
        print('{:20s}{:15.10f}'.format('SCF Energy:',self.rhfe))
        print('{:20s}{:15.10f}'.format('MP2 Correction:',self.mp2e))
        print('-'*35)
        print('{:20s}{:15.10f}'.format('FINAL ENERGY:',self.E))
        print('-'*35)
        print('='*78)


""" 
ALTERNATE INTEGRAL TRANSFORMATION ALGORITHMS
--------------------------------------------
1) ONE LINE EINSUM; 28.3423497677 seconds for H2O cc-pVDZ
gt = np.einsum('ijkl,iu,jv,kp,ls->uvps',g,c_occ,c_occ,c_vir,c_vir)

2) EIGHT FOR LOOPS; >2 hours for H2O cc-pVDZ
gt = np.zeros((nocc,nocc,nvir,nvir))
for i in range(nocc):
    for j in range(nocc):
        for a in range(nvir):
            for b in range(nvir):
                for mu in range(nbf):
                    for nu in range(nbf):
                        for rho in range(nbf):
                            for sigma in range(nbf):
                                gt[i,j,a,b] += g[mu,nu,rho,sigma]*c_occ[mu,i]*c_occ[nu,j]*c_vir[rho,a]*c_vir[sigma,b]
"""






        




        

    
