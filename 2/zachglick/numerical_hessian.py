import sys
sys.path.append('../../0/zachglick')
sys.path.append('../../1/zachglick')
import os
from subprocess import call

from molecule import Molecule
import hessian_to_freqs


def generate_inputs() :
    pass

def run_jobs() :
    pass

def build_hessian() :
    pass

if __name__ == '__main__':

    with open('../../extra-files/molecule.xyz', 'r') as f:
        molecule = Molecule(f.read())


    formatted_atoms = [atom[0] + str(atom[1][0]) + str(atom[1][1]) + str(atom[1][2]) for atom in molecule]
    print(formatted_atoms)
    for atom in molecule:
        print(atom[1][0])

    with open('../../extra-files/template.dat', 'r') as f:
        print((f.read()).format("test"))
    call(["psi4", "project2_input.dat"])
