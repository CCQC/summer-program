import numpy as np


class Molecule:
    """
    A Molecule class that implements the basic properties of a molecule
    """

    def __init__(self, geom_str, units="Angstrom"):
        """
        Make a new molecule
        :param geom_str: a string of the geometry
        :param units: a string of the units used (Bohr or Angstrom)
        """
        self.read(geom_str)
        self.units = units

    def __str__(self):
        """
        Format the molecule in a nice way
        """
        line_form = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
        out = "{:d}\n{:s}\n".format(len(self), self.units)
        for atom, xyz in zip(self.atoms, self.geom):
            out += line_form.format(atom, *xyz)

        return out

    def __len__(self):
        """
        :return: the length of the molecule
        """
        return len(self.geom)

    def read(self, geom_str):
        """
        Reads a geometry string and saves the results
        :param geom_str:
        
        Generates the following object variables
        :self.atoms: a list of the atoms
        :self.geom: a numpy array with the first index corresponding to the atom
                    and the second corresponding to the x, y, and z coordinates
        """
        geom = []
        self.atoms = []
        lines = geom_str.split("\n")
        num = int(lines[0])
        for line in lines[2:]:
            if line.strip() == "":
                continue
            atom, x, y, z = line.split()
            self.atoms.append(atom)
            geom.append([float(x), float(y), float(z)])

        if not num == len(geom):
            raise Exception("XYZ file no formatted correctly.\n" + geom_str)

        self.geom = np.array(geom)

    def to_bohr(self):
        """
        Convert the geometry to units of Bohr
        """
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.geom *= 1.889725989

    def to_angstrom(self):
        """
        Convert the geometry to units of Angstrom
        """
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.geom /= 1.889725989


if __name__ == "__main__":
    mol = Molecule(open("../../1/extra-files/molecule.xyz").read())
    print(mol)
    mol.to_bohr()
    print(mol)
