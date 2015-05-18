#!/usr/bin/env python3

import sys
sys.path.insert(0, '../../0/jevandezande')
from molecule import Molecule
import numpy as np
import subprocess
import os
from itertools import combinations_with_replacement as cwr

# Create a molecule
mol = Molecule(open('../../1/extra-files/molecule.xyz').read(), 'Bohr')

single_displace_form = 'disp/a{}_{}{}'
double_displace_form = 'disp/a{}_{}{}.a{}_{}{}'

input_form = """molecule {{
{}
units bohr
}}

set {{
    basis 3-21G
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
    Move the specified coordinate and atom the specified amount
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


def get_energy(output_file):
    lines = open(output_file).readlines()
    for line in reversed(lines):
        if line[:18] == '    Total Energy =':
            return float(line.split()[-1])

def E1(atom, coord, sign):
    output_file = single_displace_form.format(atom, coord, sign) + '.out'
    return get_energy(output_file)

def E2(atom1, coord1, atom2, coord2, sign1, sign2):
    if atom1 > atom2:
        atom1, atom2 = atom2, atom1
    if atom1 == atom2 and coord1 > coord2:
        coord1, coord2 = coord2, coord1

    output_file = double_displace_form.format(atom1, coord1, sign1,
                                              atom2, coord2, sign2) + '.out'
    
    return get_energy(output_file)


DISP = 0.005
class Calc_Iter:
    def __init__(self, mol):
        self.mol = mol
        
    def __iter__(self):
        # Iterate again over all possible unique combinations of atoms and their coordinates
        # Allows for double counting of an atom
        for (atom1, coord1), (atom2, coord2) in cwr(np.ndindex(self.mol.geom.shape), r=2):
            print(atom1, coord1, "    ", atom2, coord2)
            # Iterate over all ways to displace, don't displace the same atom coordinate twice
            if atom1 == atom2 and coord1 == coord2:
                yield atom1, coord1, atom2, coord2, DISP, 0
                yield atom1, coord1, atom2, coord2, -DISP, 0
            else:
                yield atom1, coord1, atom2, coord2, DISP, DISP
                yield atom1, coord1, atom2, coord2, -DISP, -DISP


os.makedirs('disp/', exist_ok=True)

#mol = Molecule('1\n\nH 0 0 0')

e0_in = 'disp/e0.in'
write_input(e0_in, mol)
run(e0_in)
lines = open('disp/e0.out').readlines()
for line in reversed(lines):
    if line[:18] == '    Total Energy =':
        e0 = float(line.split()[-1])
print("E0: " + str(e0))

# Iterate over all atoms and their coordinates
for atom1, coord1 in np.ndindex(mol.geom.shape):
    # Iterate again over all atoms and their coordinates
    my_iter = Calc_Iter(mol)
    for atom1, coord1, atom2, coord2, disp1, disp2 in my_iter:
        sign1 = '+' if disp1 > 0 else '-'
        if disp2 == 0:
            inp_name = single_displace_form.format(atom1, coord1, sign1) + '.in'
        else:
            sign2 = '+' if disp2 > 0 else '-'
            inp_name = double_displace_form.format(atom1, coord1, sign1, atom2, coord2, sign2) + '.in'
        displaced = displace(mol, atom1, coord1, disp1, atom2, coord2, disp2)
        write_input(inp_name, displaced)
        run(inp_name)


H = np.zeros((3*len(mol), 3*len(mol)))

for atom1, coord1 in np.ndindex(mol.geom.shape):
    for atom2, coord2 in np.ndindex(mol.geom.shape):
        if atom1 == atom2 and coord1 == coord2:
            print(atom1, coord1)
            H[3*atom1 + coord1][3*atom2 + coord2] = \
                (E1(atom1, coord1, '+') + E1(atom2, coord2, '-') - e0)/(DISP**2)
        else:
            print(atom1, coord1, atom2, coord2)
            H[3*atom1 + coord1][3*atom2 + coord2] = \
                (E2(atom1, coord1, atom2, coord2, '+', '+')
                + E2(atom1, coord1, atom2, coord2, '-', '-')
                - E1(atom1, coord1, '+') - E1(atom1, coord1, '-')
                - E1(atom2, coord2, '+') - E1(atom2, coord2, '-')
                + 2*e0) / (2*DISP**2)


hess_form = '{: 5.3f} '*len(mol)*3
out = '\n'.join([hess_form.format(*row) for row in H])
print(out)
