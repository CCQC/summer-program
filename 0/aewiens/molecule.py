#! /usr/bin/python

import numpy as np
from scipy import linalg as la

class Molecule:
    """
    implements the basic properties of a molecule
    """

    def __init__(self, geom_str, units="Angstrom"):
        self.read(geom_str)
        self.units = units

    def read(self, geom_str):
        """
        Read in the geometry (as a list of strings, by readlines() ), and store 2 class variables:
        1. self.atoms (a list of labels of the N atoms)
        2. self.geom (an Nx3 numpy array of the x,y,z coordinates of each of the atoms) 
        """
        self.atoms = [] 
        geom = []
        for line in geom_str.split('\n')[2:]:
            if line.strip() == '': 
                continue
            atom, x, y, z = line.split()[:4]
            self.atoms.append(atom)
            geom.append([float(x),float(y),float(z)])
        self.geom = np.array(geom)

    def __str__(self):
        """
        Print the molecule in a nice format
        """
        out = "{:d}\n{:s}\n".format(len(self),self.units)
        for atom, xyz in zip(self.atoms, self.geom):
            out += "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n".format(atom, *xyz)
        return out

    def __len__(self):
        """
        return the length of the molecule (think of as the # of atoms, N)
        """
        return len(self.geom)


    def bohr(self):
        """
        return the geometry in bohr
        """
        if self.units == "Angstrom":
            self.geom *= 1.889725989
            self.units == "Bohr"
        return self.geom

    def angs(self):
        """
        return the geometry in Angstroms
        """
        if self.units == "Bohr":
            self.geom /= 1.889725989
            self.units = "Angstrom"
        return self.geom

    def copy(self):
        return Molecule(str(self),self.units)

#### Example of an input
"""
f = open("../../extra-files/molecule.xyz","r")
f = f.read()
"""
f = open("TSTSTS.xyz", "r")
f = f.read()

mol = Molecule(f,"Bohr")
mol.angs()
g = open("test.xyz","w")
g.write( str(mol)) 
g.close()

#print(mol.geom)
#print(mol.atoms)
#print( mol.__len__() )
#print( mol.__str__() )
#print( mol.angs() )
#print( mol.bohr() )
#print( mol.copy() )
