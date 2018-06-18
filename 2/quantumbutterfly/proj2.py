#################################################################
# Numerical Hessian - written by quantumbutterfly 2015-07-24
# Last updated: 2015-07-28
#################################################################

import os
import sys
sys.path.append("../../0/quantumbutterfly/")
from molecule import Molecule
import numpy as np

#Initialize a molecule object from .xyz molecular geom file
geom_str = open("../../extra-files/molecule.xyz").read()
mol = Molecule(geom_str)
#mol.mol_print()
#Convert geometry to Bohr
mol = mol.ang_to_bohr()

#Copy molecular attributes
num = mol.num
atoms = mol.atoms
coords = mol.coords

#Displacement in Bohr
h = 0.005

#Function to construct input files given a directory, list of atom
#labels and a list of lists of coordinates
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
fl = coords.flat #Element-by-element iterator
for i in range(0, 3*num):
    #ith element displacement forward
    if not os.path.exists("f" + str(i) + "__"):
        os.mkdir("f" + str(i) + "__")
    f_coords = coords.flatten()
    f_coords[i] = f_coords[i] + h
    if not os.path.exists("f" + str(i) + "__/input.dat"):
        make_input("f" + str(i) + "__",atoms,np.reshape(f_coords,(num, 3)).tolist())

    #ith element displacement backward
    if not os.path.exists("b" + str(i) + "__"):
        os.mkdir("b" + str(i) + "__")
    b_coords = coords.flatten()
    b_coords[i] = b_coords[i] - h
    if not os.path.exists("b" + str(i) + "__/input.dat"):
        make_input("b" + str(i) + "__",atoms,np.reshape(b_coords,(num, 3)).tolist())

	#Double displacement directory loop
    for j in range(i + 1, 3*num):
        if not os.path.exists("f" + str(i) + "f" + str(j)):
            os.mkdir("f" + str(i) + "f" + str(j))
        ff_coords = coords.flatten()
        ff_coords[i] = ff_coords[i] + h
        ff_coords[j] = ff_coords[j] + h
        if not os.path.exists("f" + str(i) + "f" + str(j) + "/input.dat"):
            make_input("f" + str(i) + "f" + str(j),atoms,np.reshape(ff_coords,(num, 3)).tolist())

        if not os.path.exists("b" + str(i) + "b" + str(j)):
            os.mkdir("b" + str(i) + "b" + str(j))
        bb_coords = coords.flatten()
        bb_coords[i] = bb_coords[i] - h
        bb_coords[j] = bb_coords[j] - h
        if not os.path.exists("b" + str(i) + "b" + str(j) + "/input.dat"):
            make_input("b" + str(i) + "b" + str(j),atoms,np.reshape(bb_coords,(num, 3)).tolist())
