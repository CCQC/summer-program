#!/usr/bin/env python

import numpy as np
import sys

sys.path.insert(0,'../../extra-files/')
from masses import get_mass, get_charge

class Molecule:
    """
    object that contains the basic properties of a molecule
    """

    def __init__(self, geom_str, units="Angstrom"):

        self.units = units
        self.read(geom_str)
        self.masses = [ float(get_mass(i)) for i in self.atoms ]
        self.charges = [ int(get_charge(i)) for i in self.atoms ]

    def read(self, geom_str):
        """
        Read in the geometry and store 2 class variables:
        :self.atoms: a list of labels of the N atoms
        :self.geom:  an Nx3 2darray of x,y,z 
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

    def __len__(self):
        """
        :return: number of atoms
        """
        return len(self.geom)

    def __str__(self):
        """
        Print the molecule in an xyz file format
        """
        out = "{:d}\n{:s}\n".format(len(self),self.units)
        for atom, xyz in zip(self.atoms, self.geom):
            out += "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n".format(atom, *xyz)
        return out

    def bohr(self):
        """
        :return: ndarray of the Cartesian geometry in bohr
        """
        if self.units == "Angstrom":
            self.geom *= 1.889725989
            self.units = "Bohr"
        return self.geom

    def angs(self):
        """
        :return: ndarray of the Cartesian geometry in Angstroms
        """
        if self.units == "Bohr":
            self.geom /= 1.889725989
            self.units = "Angstrom"
        return self.geom

    def copy(self):
        """
        :return: a fresh copy of the molecule object
        """
        return Molecule(str(self),self.units)


if __name__ == '__main__':
    f = open("../../extra-files/molecule.xyz","r").read()
    mol = Molecule(geom_str,"Angstrom")
    mol.to_bohr()
