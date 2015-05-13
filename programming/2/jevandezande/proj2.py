#!/usr/bin/env python3

import sys
sys.path.insert(0, '../../0/jevandezande')
from molecule import Molecule
import numpy as np
import subprocess
import os

# Create a molecule
mol = Molecule(open('../../1/extra-files/molecule.xyz').read(), 'Bohr')


input_form = """molecule {{
{}
}}

set {{
    basis cc-pVDZ
}}

energy('{}')
"""


def write_input(name, molecule, theory='scf'):
    """
    Write an input file with the specified name
    :param dirname: the directory to write the input file in
    :param geom_str: the molecular geometry as a string
    :param theory: the level of theory you would like to use
    """
    # Remove the first two lines
    geom_str = '\n'.join(str(molecule).split('\n')[2:])
    # Write
    open(name, 'w').write(input_form.format(geom_str, theory))


def displace(molecule, atom1, coord1, disp1, atom2, coord2, disp2):
    """
    """
    copy = molecule.copy()
    copy.geom[atom1][coord1] += disp1
    copy.geom[atom2][coord2] += disp2
    return copy

def run(name):
    """
    Run the specified input file
    """
    # Name the oputput file the same, except for the extension
    out = '.'.join(name.split('.')[:-1]) + '.out'
    # Run
    subprocess.check_call('. ~/setup.zsh; psi4 -i {} -o {}'.format(name, out), shell=True)

os.makedirs('disp/', exist_ok=True)

name_form = 'disp/input.a{}_{}r{}.a{}_{}r{}.in'
# Iterate over all atoms and their coordinates
for atom1, coord1 in np.ndindex(mol.geom.shape):
    # Iterate again over all atoms and their coordinates
    for atom2, coord2 in np.ndindex(mol.geom.shape):
        for sign1, disp1 in [('+', 0.005), ('-', 0.005)]:
            for sign2, disp2 in [('+', 0.005), ('-', 0.005)]:
                inp_name = name_form.format(atom1, sign1, coord1, atom2, sign2, coord2)
                displaced = displace(mol, atom1, coord1, disp1, atom2, coord2, disp2)
                write_input(inp_name, mol)
                run(inp_name)
