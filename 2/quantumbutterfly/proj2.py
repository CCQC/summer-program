#################################################################
# Numerical Hessian - written by quantumbutterfly 2015-07-24
# Last updated: 2015-07-24
#################################################################

import os
import sys
sys.path.append("../../0/quantumbutterfly/")
from molecule import Molecule
import numpy as np

#Displacements fixed to 0.005 Bohr
#1 directory for the reference config
#2*3N directories for single displacements

#From the molecular geometry file, initialize a molecule object
geom_str = open("../../extra-files/molecule.xyz").read()
mol = Molecule(geom_str)
#mol.mol_print()
#Convert geometry to Bohr to match units in input
mol = mol.ang_to_bohr()

#Make a copy of the molecular coordinates for generating displacement
#geometries.
num = mol.num
atoms = mol.atoms
coords = mol.coords

#Displacement, in Bohr
#Make sure the coordinates I'm displacing are also in Bohr!
h = 0.005

#Make a directory for the reference configuration
if not os.path.exists("____"):
    os.mkdir("____")

#Loop to generate directories in which to run computations
for i in range(0, 3*num):
    if not os.path.exists("f" + str(i) + "__"):
        os.mkdir("f" + str(i) + "__")
    
	#Loop to generate double displacements
    for j in range(i, 3*num):
        if not os.path.exists("f" + str(i) + "f" + str(j)):
            os.mkdir("f" + str(i) + "f" + str(j))
        if not os.path.exists("b" + str(i) + "b" + str(j)):
            os.mkdir("b" + str(i) + "b" + str(j))



#From instructions: the third argument to the make_input() function is a list of
#coordinates. I'm making this list by converting my numpy array to a list.
#Function to print an input to a file given a set of coordinates and atom labels
#Also, note that it's very important to define this function BEFORE using it!
#Alternatively, it might be useful to include the function in a separate file,
#then include it by linking that file to the main code...

#Third argument to this function is a list (of lists of coordinates)
def make_input(dirname,labels,coords):
    with open(dirname + "/input.dat",'w') as f:
        f.write("molecule {\n")
        f.write("units bohr\n")
        for k in range(0, num):
            atom = labels[k]
            x,y,z = coords[k]
            f.write("  {:<2}{:14.10f}{:14.10f}{:14.10f}\n".format(atom,x,y,z))
        #The write function needs a string, or something that acts like one. 
        #For whatever reason, my coordinates are being truncated. The list
        #of lists of coordinates is not truncating them; not sure what the
        #deal is.
        f.write("\n}\n\n")
        f.write("set basis cc-pVDZ\n")
        f.write("energy('scf')\n")


#Combine with the directory loop?
#Make an input file for the reference configuration
if not os.path.exists("____/input.dat"):
    make_input("____",atoms,coords)

#Loop for single-displacement input files
#Still need to write code that will actually displace the coordinate
for i in range(0, 3*num):
    if not os.path.exists("f" + str(i) + "__/input.dat"):
        #Not the actual syntax for making a directory; input real code
        coords.tolist()
        make_input("f" + str(i) + "__",atoms,coords)
    
	#Loop to generate double displacement input files
	#Still need to write code to displace the coordinates
    for j in range(i, 3*num):
        if not os.path.exists("f" + str(i) + "f" + str(j) + "/input.dat"):
            make_input("f" + str(i) + "f" + str(j),atoms,coords)
        if not os.path.exists("b" + str(i) + "b" + str(j) + "/input.dat"):
            make_input("b" + str(i) + "b" + str(j),atoms,coords)
