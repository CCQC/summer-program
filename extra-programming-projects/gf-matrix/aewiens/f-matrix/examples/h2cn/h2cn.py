#!/usr/bin/env python

import sys, numpy as np

sys.path.append("../../")
from molecule import Molecule
from frequencies import Frequencies

f    = open("template.dat","r").read()
h2cn = Molecule(f)
h2cn.toBohr()
h2cn.toRadian()

##  from CFOUR, CCSD(T) / cc-pVTZ
hessianString = open("FCMINT","r").read()

##  from Andreas' 8950C code 
G = np.array([[  1.54746225e-01,  -4.31856759e-02,  -6.50791831e-02,  -4.31856759e-02,
   -6.50791831e-02,   0.00000000e+00],
 [ -4.31856759e-02,   1.07556906e+00,  -5.68642680e-02,  -3.85732708e-02,
    1.24315992e-01,   1.35385841e-17],
 [ -6.50791831e-02,  -5.68642680e-02,   1.05825984e+00,   1.24315992e-01,
   -1.29274488e-01,  -2.34948504e-17],
 [ -4.31856759e-02,  -3.85732708e-02,   1.24315992e-01,   1.07556906e+00,
   -5.68642680e-02,   1.35385841e-17],
 [ -6.50791831e-02,   1.24315992e-01,  -1.29274488e-01,  -5.68642680e-02,
    1.05825984e+00,  -2.34948504e-17],
 [  0.00000000e+00,   1.35385841e-17,  -2.34948504e-17,   1.35385841e-17,
   -2.34948504e-17,   3.13101008e+00]])


freq = Frequencies(h2cn,hessianString,G)
ref  = [940.8297,979.3300,1386.1199,1668.3542,3001.4589,3068.8091][::-1]


print("----------------------------------------")
for i,j in enumerate( freq.getFrequencies() ):
	print("{:d}. {:15.12f} {:20.12f}".format(i+1,j,ref[i]) )
print("----------------------------------------")
