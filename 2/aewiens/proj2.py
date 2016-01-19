#!/usr/bin/env python

import os
import re
import sys
sys.path.insert(0, '../../0/aewiens/')

from molecule import Molecule
import numpy as np


mol = Molecule(open("../../extra-files/molecule.xyz","r").read(),"Bohr")
N = len(mol)
h = 0.005

input_file = """molecule {{
{}
units bohr
}}

set {{
    basis cc-pVDZ
    e_convergence 12
    scf_type pk
    }}

energy('scf')
"""

def make_input(dirname,atoms,geom):
    """
    Write a psi4 input file with specified name
    :param dirname: directory where we'll write the input file
    :param atoms: list of strings of the atom labels for the molecule
    :param geom: numpy matrix of the xyz coordinates for the molecule
    """
    os.mkdir("%s" % dirname)
    xyz = [('  ').join([atoms[i]] + [str(geom[i,j]) for j in range(3)]) for i in range(N)]
    open("{:s}/input.dat".format(dirname), 'w').write(input_file.format(('\n').join(xyz)))


def run_input(dirname):
    """
    Run a psi4 single-point energy with specified name
    :param dirname: directory where we'll run the input file
    """
    os.chdir(dirname)
    os.system('psi4')
    os.chdir('..')


def E(i,j,hi,hj):
    """
    Pull energy from output.dat and throw an error if not found
    :param i: index of atom 1
    :param j: index of atom 2
    :param hi: displacement of atom 1 (-1, 0, or 1, corresponds to -h, 0, or h)
    :param hj: displacement of atom 2
    """
    dirname = "X%dX%d_%d%d" % (i,j,hi,hj)
    out_str = open("%s/output.dat" % dirname, "r").read()
    match = re.findall("Total Energy\s=\s+-\d+.\d+",out_str)
    if match == []:
        out = "Cannot find energy!"
    else:
        out = float(match[0].split()[-1])
    return out


####  Run reference configuration    ####

make_input("X0X0_00",mol.atoms,mol.geom)
run_input("X0X0_00")


####   Run single displacements   ####

for i in range(3*N):
   forward = "X%dX0_10" % i
   reverse = "X%dX0_-10" % i
   geom_copy = mol.copy().geom
   geom_copy[i/3,i%3] +=h

   make_input(forward,mol.atoms,geom_copy)

   geom_copy[i/3,i%3] -= 2*h
   make_input(reverse,mol.atoms,geom_copy)

   run_input(forward)
   run_input(reverse)


####   Run double displacements    ######

for i in range(3*N):
    for j in range(i):
        forward = "X%dX%d_11" % (i,j)
        reverse = "X%dX%d_-1-1" % (i,j)
        geom_copy2 = mol.copy().geom

        geom_copy2[i/3,i%3] += h
        geom_copy2[j/3,j%3] += h

        make_input(forward,mol.atoms,geom_copy2)

        geom_copy2[i/3,i%3] -= 2*h
        geom_copy2[j/3,j%3] -= 2*h

        make_input(reverse,mol.atoms,geom_copy2)

        run_input(forward)
        run_input(reverse)


E0 = E(0,0,0,0)
H = np.zeros((3*N,3*N))

for i in range(3*N):
    H[i,i]= (E(i,0,1,0)+E(i,0,-1,0)-2*E0)/(h**2)
    for j in range(0,i):
        H[i,j] = H[j,i] = (E(i,j,1,1)+E(i,j,-1,-1)-E(i,0,1,0)-E(j,0,1,0)-E(j,0,-1,0)-E(i,0,-1,0)+2*E0)/(2*h**2)


np.savetxt("hessian.txt",H,"%15.7f"," ","\n")
