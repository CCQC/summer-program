#################################################################
# Numerical Hessian - written by quantumbutterfly 2015-07-24
# Last updated: 2015-07-24
#################################################################

import os
import sys
sys.path.append("../../0/quantumbutterfly/")
import molecule from molecule
import numpy as np

#Displacements fixed to 0.005 Bohr
#1 directory for the reference config
#2*3N directories for single displacements
#

#From the molecular geometry file, initialize a molecule object
geom_str = open("../../extra-files/molecule.xyz").read()
mol = Molecule(geom_str)
#Convert geometry to Bohr to match units in input
mol = mol.ang_to_bohr()
mol.print()

#Make a copy of the molecular coordinates for generating displacement
#geometries.
num = mol.num
atoms = mol.atoms
coords = mol.coords


#Loop to generate directories for running computations in
#We have 3N coordinates, but ordered starting with zero, so
#iterate from 0 to 3N-1.
for i in range(0, 3*num):
    if not path.exists("comps/f" + i + "__"):
        #Not the actual syntax for making a directory; input real code
        mkdir "f" + i + "__"

    
    for j in range(i, 3*num):
        if not path.exists("comps/f" + i + "f" + j):
            mkdir "f" + i + "f" + j
        if not path.exists("comps/f" + i + "b" + j):
            mkdir "f" + i + "b" + j
        if not path.exists("comps/b" + i + "f" + j):
            mkdir "b" + i + "f" + j
        if not path.exists("comps/b" + i + "b" + j):
            mkdir "b" + i + "b" + j




#Each time a geometry is generated in the loop structure, make sure the geometry
#is stored as the array geom, then pass it (and the set of atom labels from the
#beginning of initialization, above) to the make_input function to make an input.




#Function to print an input to a file given a set of coordinates and atom labels
def make_input(dirname,labels,coords):
    f = open("input.dat","w")
    open(input.dat).write()
    print("molecule {\n")
    print("units bohr\n")
    print(first symbol  coords\n)
    print(other symbol  coords\n)
    print("}\n\n")
    print("set basis cc-pVDZ")
    print("energy('scf')")
