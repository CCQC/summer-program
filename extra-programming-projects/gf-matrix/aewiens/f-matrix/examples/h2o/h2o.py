#!/usr/bin/env python

import sys, numpy as np

sys.path.insert(0,'../../')
from molecule import Molecule
from hessian import Hessian
from frequencies import Frequencies

ref = [1853.0679, 2335.9382, 2475.0096][::-1]

f    = open("template.dat","r").read()
mol  = Molecule(f)
mol.toBohr()
mol.toRadian()
water = Hessian(mol,f)

##  build Gmatrix for water  ## 
muO  = 1/mol.masses[0]
muH  = 1/mol.masses[1]
r    = mol.coords[0]
phi  = mol.coords[2]

G = np.zeros((3,3))
G[0,0] = G[1,1] = muO + muH
G[1,0] = G[0,1] = muO*np.cos(phi)
G[0,2] = G[2,0] = -muO*np.sin(phi)/r
G[2,1] = G[1,2] = -muO*np.sin(phi)/r
G[2,2] = 2*(muO + muH - muO*np.cos(phi))/r**2


h    = open("FCMINT","r").read()
freq = Frequencies(mol,h,G)


print("----------------------------------------")
print("    H2O Harmonics:         Reference:")
for i,j in enumerate(freq.getFrequencies()):
	print("{:d}. {:15.12f} {:20.12f}".format(i+1,j,ref[i]) )
print("----------------------------------------")
