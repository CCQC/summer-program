#!/usr/bin/python

import os
import numpy as np

def generate_inputs(mol, template, disp_size = 0.05, directory = "DISPS"):
    
    mol = mol.copy()
    mol.to_bohr()

    def write_input(coord1, coord2, k1, k2, disp_mol):
        path = directory + "/x{:d}x{:d}_{:d}{:d}".format(coord1, coord2, k1, k2)
        if not os.path.isdir(path): os.makedirs(path)
        inpf = open(path + "/input.dat",'w+')
        inpf.write(template.format(str(disp_mol)))
        inpf.close()

if __name__ == "__main__":
    import sys


