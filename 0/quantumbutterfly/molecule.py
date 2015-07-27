import numpy as np

class Molecule:
    """
    A molecule class that implements labeled geometry and number of atoms
    """

    def __init__(self, geom_str, units="Angstrom"):
        """
        Make a new molecule
        :param geom_str: a string of the desired geometry
        :param units: units of the geometry, Bohr or Angstrom
        """
        self.units = units
        self.read(geom_str)

    def read(self, geom_str):
        """
        Read in the molecular geometry
        :param geom_str: string containing molecular geometry
        """
        lines = geom_str.strip().split("\n")
        num = int(lines[0])
        atoms = []
        coords = []

        for line in lines[2:]:
            atom, x, y, z = line.split()
            atoms.append(atom)
            coords.append([float(x), float(y), float(z)])

        self.num = num
        self.atoms = atoms
        self.coords = np.array(coords)
        """
    def to_bohr(self):
       #
       # Convert geometry to Angstrom if in Bohr
       #
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.coords *= 1.889725989

    def to_ang(self):
	#   Convert geometry to Bohr if in Angstrom
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.coords /= 1.889725989
"""
    def print(self):
        print(self.num)
        print(self.atoms)
        print(self.coords)

if __name__ == "__main__":
    geom_str = open("../../extra-files/molecule.xyz").read()
    mol = Molecule(geom_str)
    mol.print()
