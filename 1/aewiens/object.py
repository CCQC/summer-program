#!/usr/bin/env python

import numpy as np
import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/aewiens')

import masses
from molecule import molecule

class Frequencies:

    def __init__(self,mol,hessian):

        self.N = mol.__len__()
        m = []
        for i in mol.atoms:
            m += [1/masses.get_mass(i))**0.5]*3
        M = np.diag(m)

        
        self.mH = M*hessian*M
        self.e, self.l = np.linalg.eigh(mH)

        Q = np.matrix(M)*np.matrix(l)
