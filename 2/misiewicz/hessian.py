#coding: UTF-8

import functools
import glob
import itertools
import numpy as np
import os
import re
import subprocess
import sys
sys.path.insert(0, '../../0/misiewicz')
from molecule import Molecule
sys.path.insert(1, '../../1/misiewicz')
from frequency_finder import frequencies

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):

    # Make our base directory. If it already exists, we have a problem.
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        raise Exception("Output directory {} already exists.".format(directory))

    geom = mol.geom
    labels = mol.labels

    # For each unordered pair of elements, we have to combine them. Start...
    for a, b in itertools.combinations_with_replacement(range(geom.size), 2):
        # FACT: a <= b.
        a, b = int(a), int(b)
        # Make a new directory, corresponding to the coupled pair.
        os.makedirs(os.path.join(directory, "{} {}".format(a, b)))
        # Give the coupled positive and positive movements.
        d1, d2 = displace(geom.copy(), {a, b}, disp_size)
        for i, dn in enumerate((d1, d2), start=1):
            # Write the file.
            with open(os.path.join(directory, str(a) + " " + str(b), str(i) + ".dat"), "w") as f:
                 xyz_str = format_for_xyz(labels, dn)
                 f.write(template.format(xyz_str))

    # Now add a file for the undisturbed geometry.
    os.makedirs(os.path.join(directory, "eq"))
    with open(os.path.join(directory, "eq", "0.dat"), "w") as f:
         xyz_str = format_for_xyz(labels, geom)
         f.write(template.format(xyz_str))

def to_rc(no):
    '''Convert the number of a Hessian element to (row, column).'''
    return (no // 3, no % 3)

def displace(matrix, cell_nos, displacement):
    '''Given a matrix, cell numbers, and a displacement, return the
    matrices displaced positively and negatively at those coordinates.'''
    mat_plus, mat_neg = matrix, matrix.copy()
    for cell_no in cell_nos:
        mat_plus[to_rc(cell_no)] += displacement
        mat_neg[to_rc(cell_no)] -= displacement
    return mat_plus, mat_neg

def format_for_xyz(labels, matr):
     xyz_string = "units bohr\n"
     format_str = "{:3}" + " {:15.10f}" * 3 + "\n"

     for label, row in zip(labels, matr):
         xyz_string += format_str.format(label, *row)
     return xyz_string        

def run_jobs(mol, command = "psi4", directory = "DISPS"):
    for filename in glob.iglob("{}/*/*.dat".format(directory)):
        subprocess.check_call(command + " -i '" + filename + "' -o '" + filename.replace(".dat", ".out'"), shell=True)

def get_energy(directory, anchor, sign, *args):
    '''Open the given filename and return the energy, as identified by an anchor phrase.'''
    # TODO: Store these values so there is less regex.
    folder = " ".join(str(arg) for arg in args)
    name = {"+": "1", "-": "2", "0": "0"}[sign]
    filename = os.path.join(directory, folder, name + ".out")
    with open(filename, "r") as f:
        text = f.read()
    match_strings = re.findall("{}\s*(-?\d*\.\d*)".format(anchor), text)
    if len(match_strings) != 1:
        raise Exception("Did not match exactly once in {}.".format(filename))
    return float(match_strings[0])

def build_hessian(mol, energy_prefix, disp_size = 0.005, directory = "DISPS"):
    hessian = np.empty([3*mol.natom, 3*mol.natom])
    e = functools.partial(get_energy, directory, energy_prefix)
    for a, b in itertools.combinations_with_replacement(range(3*mol.natom), 2):
        # A diagonal element.
        if a == b:
            hessian[a, b] = 1 / disp_size ** 2 * (
                e("+", a, a) + e("-", a, a) - 2*e("0", "eq"))
        # An off-diagonal element.
        else:
            hessian[a, b] = hessian[b, a] = 1/(2 * disp_size ** 2) * (
                e("+", a, b) + e("-", a, b) - e("+", a, a) - e("-", a, a) - e("+", b, b) - e("-", b, b) + 2 * e("0", "eq"))
    return hessian

if __name__ == "__main__":
    water = Molecule("../../extra-files/molecule.xyz")
    water.to_bohr()
    # You can put your displacement in whatever unit you like, but it has to
    # match the units of your energy. Psi4 sends energy in Hartrees, so Bohr
    # units are necessary for our displacement.
    with open("../../extra-files/template.dat", "r") as f:
        template = f.read()
    generate_inputs(water, template)
    run_jobs(water)
    hessian = build_hessian(water, "@DF-RHF Final Energy:")
    np.savetxt("hessian.out", hessian)
    frequencies(water, hessian)
