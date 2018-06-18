import os
import numpy as np
import re

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):

    mol = mol.copy()
    mol.ang_to_bohr() 

    
