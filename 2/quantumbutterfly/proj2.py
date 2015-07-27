#################################################################
# Numerical Hessian - written by quantumbutterfly 2015-07-24
# Last updated: 2015-07-27
#################################################################

import os
import sys
sys.path.append("../../0/quantumbutterfly/")
from molecule import Molecule
import numpy as np

#From the molecular geometry file, initialize a molecule object
geom_str = open("../../extra-files/molecule.xyz").read()
mol = Molecule(geom_str)
#mol.mol_print()
#Convert geometry to Bohr
mol = mol.ang_to_bohr()

#Copy molecular coordinates for displacement geometries.
num = mol.num
atoms = mol.atoms
coords = mol.coords

#Displacement, in Bohr
h = 0.005

#Function to construct input files given a directory, list of atom
#labels and a list of coordinates
def make_input(dirname,labels,coords):
    with open(dirname + "/input.dat",'w') as f:
        f.write("molecule {\nunits bohr\n")
        for k in range(0, num):
            atom = labels[k]
            x,y,z = coords[k]
            f.write("  {:<2}{:14.10f}{:14.10f}{:14.10f}\n".format(atom,x,y,z))
        f.write("\n}\n\nset basis cc-pVDZ\nenergy('scf')\n")

#Reference configuration
if not os.path.exists("____"):
    os.mkdir("____")
if not os.path.exists("____/input.dat"):
    make_input("____",atoms,coords)

#Single displacement directory loop
for i in range(0, 3*num):
    if not os.path.exists("f" + str(i) + "__"):
        os.mkdir("f" + str(i) + "__")
    if not os.path.exists("f" + str(i) + "__/input.dat"):
        coords.tolist()
        make_input("f" + str(i) + "__",atoms,coords)

	#Double displacement directory loop
    for j in range(i, 3*num):
        if not os.path.exists("f" + str(i) + "f" + str(j)):
            os.mkdir("f" + str(i) + "f" + str(j))
        if not os.path.exists("f" + str(i) + "f" + str(j) + "/input.dat"):
            make_input("f" + str(i) + "f" + str(j),atoms,coords)

        if not os.path.exists("b" + str(i) + "b" + str(j)):
            os.mkdir("b" + str(i) + "b" + str(j))
        if not os.path.exists("b" + str(i) + "b" + str(j) + "/input.dat"):
            make_input("b" + str(i) + "b" + str(j),atoms,coords)
