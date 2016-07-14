import numpy as np
import sys
import os

sys.path.insert(0, '../../0/oliviabern')
from molecule import Molecule
sys.path.insert(0, '../../extra-files')
from masses import mass

#build molecule object
geom_string = open('../../extra-files/molecule.xyz').read()
mol = Molecule(geom_string)

#build input file template
template = open('../../extra-files/template.dat').read()

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):
    m = str(mol)
    N = int(m[0])
    p = '/Users/oliviambernstein/Desktop/ccqc/'
    ph = p + 'hessian/'
    if not os.path.exists(ph):
        os.mkdir(ph)
    if not os.path.exists(ph + 'geo/'):
        os.mkdir(ph + 'geo/')
        file = open(ph + 'geo/input.dat','w')
        file.write('hello')
        file.close()


generate_inputs(mol, template)

