__author__ = 'adabbott'

import numpy as np

class Molecule(object):
    """
    A simple Molecule
    """
    def __init__(self, geom_str, units="Angstrom"):
        """
        Make a new molecule
        :param geom_str: a string of the desired geometry
        :param units: the units of the geometry, Bohr or Angstrom
        """
        self.units = units
        self.read(geom_str)

    def read(self, geom_str):
        """
        Read the molecular geometry
        :param geom_str: A string containing molecular geometry
        """
        lines = geom_str.strip().split("\n")
        num = int(lines[0])

        atoms = []
        coords = []

        for line in lines[2:]:
            atom, x, y, z = line.split()
            atoms.append(atom)
            coords.append([float(x),float(y),float(z)])
        self.atoms = atoms
        self.coords = np.array(coords)
        print("Atoms Present: ")
        print(atoms)
        print("Atomic Coordinates:")
        print(self.coords)

    def __len__(self):
        """
        Gives the number of atoms in the molecule
        :return:
        """
        return len(self.atoms) #redefines length function


    def ang_to_bohr(self):
        """
         Converts units from Angstrom to Bohr
        """
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.coords *= 1.889725989
            print("These are the units in Bohr: ")


    def bohr_to_ang(self):
        """
        Converts bohr to ang
        """
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.coords /= 1.889725989
            print("These are the units in Angstroms: ")

if __name__ == "__main__":
    mol = Molecule(open("../../extra-files/molecule.xyz").read())
    print("The length of this molecule is" ,len(mol))
    print(mol.ang_to_bohr())
    print(mol.bohr_to_ang())



