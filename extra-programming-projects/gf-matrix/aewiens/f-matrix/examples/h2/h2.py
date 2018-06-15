#!/usr/bin/env python

import sys, numpy as np

sys.path.insert(0,'../../')
from molecule import Molecule
from hessian import Hessian
from frequencies import Frequencies

f    = open("template.dat","r").read()
mol  = Molecule(f)
mol.toBohr()
mol.toRadian()

hess  = Hessian(mol,f)
hess.writeHessian()

##  build Gmatrix for H2  ## 

muH  = 1/mol.masses[0]

G      = np.matrix(np.zeros((1,1)))
G[0,0] = 2*muH

h    = open("hessian.dat","r").read()
freq = Frequencies(mol,h,G)

for i in freq.getFrequencies():
	print(i)
