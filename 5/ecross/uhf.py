#!/anaconda/bin/python

import psi4 
import numpy as np
import re
from scipy import linalg as la

class UHF(object):
    def __init__(self, mol, mints):
        self.mol, self.mints = mol, mints
        self.nbf             = self.mints.basisset().nbf()
        self.maxiter         = psi4.get_global_option('MAXITER')
        self.e_convergence   = psi4.get_global_option('E_CONVERGENCE')
        self.Vnu      = mol.nuclear_repulsion_energy()

        # Import one and two electron integrals in spatial AO basis
        self.s_spat   = np.array(mints.ao_overlap())
        self.t_spat   = np.array(mints.ao_kinetic())
        self.v_spat   = np.array(mints.ao_potential())
        self.g_spat   = np.array(mints.ao_eri())
        print(self.g_spat[:,:,:,:])

        # Transform 1- and 2-e integrals to spin AO basis
        self.s, self.t, self.v, self.g = self.spin_integrals()

        # Form orthogonalizer, X = S^(-0.5)
        x = la.sqrtm(la.inv(self.s))

        # Initial HF guess with D = 0
        self.d = np.zeros((2*self.nbf,2*self.nbf))

        # Iterate to self-consistency
        self.iterate()

        self.check_answer()

    def spin_integrals(self):
        s_spat, t_spat = self.s_spat, self.t_spat
        v_spat, g_spat = self.v_spat, self.g_spat
        nbf = self.nbf

        s = np.zeros((2*nbf,2*nbf))
        t = np.zeros((2*nbf,2*nbf))
        v = np.zeros((2*nbf,2*nbf))
        g = np.zeros((2*nbf,2*nbf,2*nbf,2*nbf))

        s[:nbf,:nbf] = s_spat[:nbf,:nbf] 
        s[nbf:,nbf:] = s_spat[:nbf,:nbf]
        t[:nbf,:nbf] = t_spat[:nbf,:nbf] 
        t[:nbf,:nbf] = t_spat[:nbf,:nbf] 
        v[nbf:,nbf:] = v_spat[:nbf,:nbf]
        v[nbf:,nbf:] = v_spat[:nbf,:nbf]
        return(s,t,v,g)

    def iterate(self):
        s, t, v, g = self.s, self.t, self.v, self.g
        

    def check_answer(self):
        output = open('output.dat','r').read()
        psi = float(re.findall('\s+Total Energy\s+\=\s+(\-\d+\.\d+)',output)[0])
        print(psi)
        
