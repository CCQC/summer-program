__author__ = 'Adam'


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
        print(atoms)
        print(coords)





if __name__ == "__main__":
    geom_str = open("../../extra-files/molecule.xyz").read()
    print(geom_str)
    mol = Molecule(geom_str)
