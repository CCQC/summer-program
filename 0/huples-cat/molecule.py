import numpy as np
import sys
sys.path.insert(0,'../../extra-files/')
from masses import get_mass, get_charge

class Molecule(object):
    def __init__(self, geom_path):
        f = open(geom_path)
        lines = f.readlines()
        f.close()
        self.natom = int(lines[0])
        self.units = lines[1]
        # initialize stuff
        self.labels = []
        self.charges = []
        self.masses = []
        self.xyz = []
        # get geoms
        for line in lines[2:]:
            atom, x, y, z = line.split()
            self.labels.append(atom)
            self.charges.append(get_charge(atom))
            self.masses.append(get_mass(atom))
            self.xyz.append((float(x), float(y), float(z)))
        self.geom = np.array(self.xyz)                            # make numpy array

    def __str__(self):
        """
        Convert into an xyz formatted string
        """
        my_str = "%i \n" % self.natom
        my_str += "%s" % self.units
        for i in range(0, self.natom):
            my_str += "%s\t%f\t%f\t%f\n" % (self.labels[i], self.xyz[i][0], self.xyz[i][1], self.xyz[i][2])
        # removes newline
        my_str = my_str[:-1]
        return my_str

    def to_bohr(self):
        if(self.units.lower() != 'bohr'):
            self.units = 'Bohr'
            self.geom *= 0.529177

    def to_angstrom(self):
        if(self.units.lower() != 'angstrom'):
            self.units = 'Angstrom'
            self.geom /= 0.529177

    def copy(self):
        return Molecule(str(self), self.units)

if __name__ == "__main__":
    mol = Molecule('../../extra-files/molecule.xyz')
    print(mol)
