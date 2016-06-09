import numpy as np

class Molecule(object):
    """
    A simple molecule
    """

    def __init__(self, geom_str, units = "Angstrom"):
        self.read(geom_string)
        self.units = units

    def read(self, geom_string):
        atoms = []
        geom = []
        lines = geom_str.split("\n")
        for item in lines[2:]:
            atom, x, y, z = item.split()
            atoms.append(atom)
            geom.append([float(x),float(y),float(z)]) 
        self.atoms = atoms
        self.geom = np.array(geom)

    def __len__(self):
        return len(atoms)

    def to_bohr(self):
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.geom *= 1.889725989

    def to_angstrom(self):
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.geom /= 1.889725989

    
    
