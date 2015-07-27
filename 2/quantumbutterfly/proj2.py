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
#

#From the molecular geometry file, initialize a molecule object
geom_str = open("../../extra-files/molecule.xyz").read()
mol = Molecule(geom_str)
mol.print()
#Convert geometry to Bohr to match units in input
#mol = mol.ang_to_bohr()
#mol.print()

#Make a copy of the molecular coordinates for generating displacement
#geometries.
num = mol.num
atoms = mol.atoms
coords = mol.coords

#Make a directory for the reference configuration
if not os.path.exists("____"):
    os.mkdir("____")
    make_input("____",atoms,coords)

#Loop to generate directories in which to run computations
for i in range(0, 3*num):
    if not os.path.exists("f" + str(i) + "__"):
        #Not the actual syntax for making a directory; input real code
        os.mkdir("f" + str(i) + "__")
    
	#Loop to generate double displacements
    for j in range(i, 3*num):
        if not os.path.exists("f" + str(i) + "f" + str(j)):
            os.mkdir("f" + str(i) + "f" + str(j))
        if not os.path.exists("f" + str(i) + "b" + str(j)):
            os.mkdir("f" + str(i) + "b" + str(j))
        if not os.path.exists("b" + str(i) + "f" + str(j)):
            os.mkdir("b" + str(i) + "f" + str(j))
        if not os.path.exists("b" + str(i) + "b" + str(j)):
            os.mkdir("b" + str(i) + "b" + str(j))


"""
#From instructions: the third argument to the make_input() function is a list of
#coordinates. I'm making this list by converting my numpy array to a list.
#Function to print an input to a file given a set of coordinates and atom labels
def make_input(dirname,labels,coords):
    f = os.open(dirname + "/input.dat",os.O_RDWR|os.O_CREAT)
    os.write(f,"molecule \n")
    os.write(f,"units bohr\n")
    os.write(f,"AT xcoord ycoord zcoord\n")
    os.write(atoms)
    os.write(coords)
    os.write(f,"}\n\n")
    os.write(f,"set basis cc-pVDZ\n")
    os.write(f,"energy('scf')\n")


#To avoid the problem of having no input file if the directory already
#exists, this input file loop is independent. Search for a way to
#combine the two loops.
#Make an input file for the reference configuration
if not os.path.exists("____/input.dat"):
    make_input("____",atoms,coords)

#Loop for single-displacement input files
#Still need to write code that will actually displace the coordinate
for i in range(0, 3*num):
    if not os.path.exists("f" + str(i) + "__/input.dat"):
        #Not the actual syntax for making a directory; input real code
        make_input("f" + str(i) + "__",atoms,coords)
    
	#Loop to generate double displacement input files
	#Still need to write code to displace the coordinates
    for j in range(i, 3*num):
        if not os.path.exists("f" + str(i) + "f" + str(j) + "/input.dat"):
            make_input("f" + str(i) + "f" + str(j),atoms,coords)
        if not os.path.exists("f" + str(i) + "b" + str(j) + "/input.dat"):
            make_input("f" + str(i) + "b" + str(j),atoms,coords)
        if not os.path.exists("b" + str(i) + "f" + str(j) + "/input.dat"):
            make_input("b" + str(i) + "f" + str(j),atoms,coords)
        if not os.path.exists("b" + str(i) + "b" + str(j) + "/input.dat"):
            make_input("b" + str(i) + "b" + str(j),atoms,coords)


#Each time a geometry is generated in the loop structure, make sure the geometry
#is stored as the array geom, then pass it (and the set of atom labels from the
#beginning of initialization, above) to the make_input function to make an input.
"""


"""
#From instructions: the third argument to the make_input() function is a list of
#coordinates. I'm making this list by converting my numpy array to a list.
#Function to print an input to a file given a set of coordinates and atom labels
def make_input(dirname,labels,coords):
    f = os.open(dirname + "/input.dat",os.O_RDWR|os.O_CREAT)
    os.write("molecule {\n")
    os.write("units bohr\n")
    os.write("AT xcoord ycoord zcoord\n")
    os.write(atoms)
    os.write(coords)
    os.write("}\n\n")
    os.write("set basis cc-pVDZ\n")
    os.write("energy('scf')\n")
"""
