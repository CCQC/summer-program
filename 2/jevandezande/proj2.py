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

single_displace_form = 'disp/q{}{}'
double_displace_form = 'disp/q{}{}_q{}{}'

input_form = """molecule {{
{}
units bohr
}}

set {{
    basis sto-3g
    e_convergence   10
    d_convergence   10
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


def displace(molecule, q1, disp1, q2, disp2):
    """
    Move the specified coordinate and atom the specified amount
    Must convert from q to xyz
    :param molecule: the molecule to be displaced
    :param q1: coordinate 1
    :param disp1: displacement for coordinate 1
    :param q2: coordinate 2
    :param disp2: displacement for coordinate 2
    :return: a new Molecule with the specified displacements
    """
    copy = molecule.copy()
    atom1 = int(q1 / 3)
    coord1 = q1 % 3
    atom2 = int(q2 / 3)
    coord2 = q2 % 3
    copy.geom[atom1][coord1] += disp1
    copy.geom[atom2][coord2] += disp2

    return copy


def run(name):
    """
    Run the specified input file
    :param name: the name of the input file
    """
    # Name the oputput file the same, except for the extension
    out = '.'.join(name.split('.')[:-1]) + '.out'
    # Run
    subprocess.check_call('psi4 -i {} -o {}'.format(name, out), shell=True)


def get_energy(output_file):
    """
    Read the energy from a file, throw an exception if it cannot find the energy
    :param output_file: the name of the file containing the energy
    :return: the energy as a float
    """
    lines = open(output_file).readlines()
    for line in reversed(lines):
        if line[:18] == '    Total Energy =':
            return float(line.split()[-1])
    raise Exception('Cannot find the energy')


def E1(q1, sign):
    """
    Read a singly displaced energy
    :param q1: the coordinate
    :param sign: the sign of the displacement
    :return: energy as a float
    """
    output_file = single_displace_form.format(q1, sign) + '.out'
    return get_energy(output_file)


def E2(q1, sign1, q2, sign2):
    """
    Read a doubly displaced energy
    :param q1: coordinate 1
    :param sign1: the sign of displacement 1
    :param q2: coordinate 2
    :param sign2: the sign of displacement 2
    :return: energy as a float
    """
    # Use symmetry
    if q1 > q2:
        q1, q2 = q2, q1

    output_file = double_displace_form.format(q1, sign1, q2, sign2) + '.out'
    
    return get_energy(output_file)


DISP = 0.005
os.makedirs('disp/', exist_ok=True)


# Run the reference geometry
e0_in = 'disp/e0.in'
write_input(e0_in, mol)
run(e0_in)
e0 = get_energy('disp/e0.out')
print("E0: " + str(e0))

# Iterate over all pairs of coordinates in the upper triangle
for q1, q2 in cwr(range(3*len(mol)), r=2):
    # Displace forwards and back
    for disp in [DISP, -DISP]:
        sign = '+' if disp > 0 else '-'
        if q1 == q2:
            print(q1, sign)
            inp_name = single_displace_form.format(q1, sign) + '.in'
            displaced = displace(mol, q1, disp, q2, 0)
        else:
            print(q1, q2, sign)
            inp_name = double_displace_form.format(q1, sign, q2, sign) + '.in'
            displaced = displace(mol, q1, disp, q2, disp)
        write_input(inp_name, displaced)
        run(inp_name)


# Form an empty Hessian matrix
H = np.zeros((3*len(mol), 3*len(mol)))

# Iterate over all coordinates in the hessian
for q1, q2 in np.ndindex(H.shape):
    # Diagonal elements
    if q1 == q2:
        H[q1][q2] = (E1(q1, '+') + E1(q1, '-') - 2*e0)/(DISP**2)
    # Off-diagonal elements
    else:
        H[q1][q2] = (E2(q1, '+', q2, '+') + E2(q1, '-', q2, '-')
            - E1(q1, '+') - E1(q1, '-') - E1(q2, '+') - E1(q2, '-')
            + 2*e0) / (2*DISP**2)


hess_form = '{: 15.13f} '*len(mol)*3
out = '\n'.join([hess_form.format(*row) for row in H])
print(out)

open('hessian.dat', 'w').write(out)
