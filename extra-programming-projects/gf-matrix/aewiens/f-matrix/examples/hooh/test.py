#!/usr/bin/env python

import sys, numpy as np

sys.path.insert(0,'../../')

from molecule import Molecule
from frequencies import Frequencies

f = open("template.dat","r").read()
mol = Molecule(f)

G = np.array(
      [[ 1.05475559, -0.01480596, -0.05945093,  0.00932732, -0.        ,         0.        ],
       [-0.01480596,  0.12503974, -0.09056411,  0.        , -0.        ,        -0.        ],
       [-0.05945093, -0.09056411,  1.24250391, -0.07336291,  0.        ,         0.        ],
       [ 0.00932732,  0.        , -0.07336291,  2.36324802,  0.        ,         0.        ],
       [-0.        , -0.        ,  0.        ,  0.        ,  1.05475559,        -0.02536097],
       [ 0.        , -0.        ,  0.        ,  0.        , -0.02536097,         1.18289273]])


h = open('hessian.dat','r').read()
freq = Frequencies(mol,h,G)

for i in freq.getFrequencies():
	print(i)

