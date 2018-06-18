#!/usr/bin/python

import sys,numpy as np
sys.path.insert(0,'/Users/avery/git/summer-program/0/aewiens/')
sys.path.insert(0,'/Users/avery/git/summer-program/1/aewiens/')
sys.path.insert(0,'/Users/avery/git/summer-program/2/aewiens/')
sys.path.insert(0,'/Users/avery/git/summer-program/extra-files/')

from molecule import Molecule
from hessian import Hessian
from frequencies import Frequencies

f  = open("input.dat","r").read()
h2o = Molecule(f)
h2o.bohr() 

hessian = Hessian(h2o,"template.dat")
hessian.write_Hessian() 


hessian = open("hessian.dat","r").read()
freq    = Frequencies(h2o,hessian)
freq.frequency_output("modes.xyz")
