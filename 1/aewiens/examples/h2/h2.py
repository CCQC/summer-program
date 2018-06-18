#!/usr/bin/python

import sys
sys.path.append("/Users/avery/git/summer-program/extra-files/")
sys.path.insert(0,'/Users/avery/git/summer-program/0/aewiens/')
sys.path.insert(0,'/Users/avery/git/summer-program/1/aewiens/')
sys.path.insert(0,'/Users/avery/git/summer-program/2/aewiens/')
from molecule import Molecule
from hessian import Hessian
from frequencies import Frequencies

f  = open("input.dat","r").read()
h2 = Molecule(f)
h2.bohr() 

hessian = Hessian(h2,"template.dat")
#hessian.write_Hessian() 

hessian = open("hessian.dat","r").read()
freq    = Frequencies(h2,hessian)
freq.frequency_output("modes.xyz")
