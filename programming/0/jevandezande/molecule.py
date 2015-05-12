import numpy as np


class Molecule:
    """
    A Molecule class that implements the basic properties of a molecule
    """

    def __init__(self, geom_str, units='Angstrom'):
        self.num, self.atoms, self.geom = Molecule.read(geom_str)
        self.units = units

    def __str__(self):
        """
        This formats the molecule in a nice way
        """
        line_form = '{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n'
        out = '{:d}\n{:s}\n'.format(self.num, self.units)
        for i in range(self.num):
            out += line_form.format(self.atoms[i], *self.geom[i])

        return out

    @staticmethod
    def read(geom_str):
        """
        Read a geometry string
        :param geom_str:
        :return: atoms_list, np.array(xyz_coords)
        """
        geom = []
        atoms = []
        lines = geom_str.split('\n')
        num = int(lines[0])
        for line in lines[2:]:
            if line.strip() == '':
                continue
            atom, x, y, z = line.split()
            atoms.append(atom)
            geom.append([float(x), float(y), float(z)])

        if not num == len(geom):
            raise Exception("XYZ file no formatted correctly.\n" + geom_str)

        return num, atoms, np.array(geom)

    def to_bohr(self):
        """
        Convert the geometry to units of Bohr
        """
        if self.units == 'Angstrom':
            self.units = 'Bohr'
            self.geom *= 1.889725989

    def to_angstrom(self):
        """
        Convert the geometry to units of Angstrom
        """
        if self.units == 'Bohr':
            self.units = 'Angstrom'
            self.geom /= 1.889725989


if __name__ == '__main__':
    mol = Molecule(open('../extra-files/molecule.xyz').read())
    print(mol)
    mol.to_bohr()
    print(mol)
