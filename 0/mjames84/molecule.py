import numpy as np
import sys
sys.path.insert(0,'../../extra-files')
from masses import get_charge, get_mass

class Molecule(object):
    """
    A simple molecule
    """

    def __init__(self, geom_str, units = "Angstrom"):
        self.read(geom_str)
        self.units = units
        self.charge = [get_charge(atom) for atom in self.atoms]
        self.mass = [get_mass(atom) for atom in self.atoms]

    def read(self, geom_str):
        atoms = []
        geom = []
        lines = geom_str.strip().split("\n")
        for item in lines[2:]:
            atom, x, y, z = item.split()
            atoms.append(atom.upper())
            geom.append([float(x),float(y),float(z)]) 
        self.atoms = atoms
        self.geom = np.array(geom)

    def __len__(self):
        return len(self.atoms)

    def to_bohr(self):
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.geom *= 1.889725989

    def to_angstrom(self):
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.geom /= 1.889725989

    def __str__(self):
        line_form = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
        out = "{:d}\n{:s}\n".format(len(self), self.units)
        for atom, xyz in zip(self.atoms, self.geom):
            out += line_form.format(atom, *xyz)
        return out 

    def copy(self):
        return Molecule(self.geom, self.units)

if __name__ == '__main__':
    geom_str = open("../../extra-files/molecule.xyz","r").read()
    mol = Molecule(geom_str)
    mol.to_bohr()
    mol.to_angstrom()
 print(mol)

    
    
