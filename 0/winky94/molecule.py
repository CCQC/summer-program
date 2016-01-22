import numpy as np

class Molecule(object):
    """
    A simple molecule
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
        Reads in the molecular geometry
        :param geom_str: a string containing molecular geometry
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

    def __len__(self):
        """

        :return: the length of the molecule
        """
        return len(self.atoms)

    def ang_to_bohr(self):
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.coords *= 1.889725989

    def bohr_to_ang(self):
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.coords /= 1.889725989

    def __str__(self):
        """
        Format the molecule in a nice way
        """
        line_form = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
        out = "{:d}\n{:s}\n".format(len(self), self.units)
        for atom, xyz in zip(self.atoms,self.coords):
            out += line_form.format(atom, *xyz)
        return out
"""
    def print(self):
        print(self.num)
        print(self.atoms)
        print(self.coords)
"""
if __name__ == "__main__":
    geom_str = open("../../extra-files/molecule.xyz").read()
    mol = Molecule(geom_str)
    mol.ang_to_bohr()
    print(mol)



