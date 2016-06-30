import numpy as np
import sys
sys.path.insert(0,'../../extra-files')
from masses import get_charge, get_mass

class Molecule(object):
    """
    A simple molecule
    """

    def __init__(self, geom, units = "Angstrom"):
        self.read(geom)
        self.units = units
        self.charge = [get_charge(atom) for atom in self.atoms]
        self.mass = [get_mass(atom) for atom in self.atoms]

    def read(self, geom):
        atoms = []
        xyz = []
        geom = open('../../extra-files/molecule.xyz','r')
        for line in geom.readlines()[2:]:
            atom, x, y, z = line.split()
            atoms.append(atom.upper())
            xyz.append([float(x),float(y),float(z)]) 
        self.atoms = atoms
        self.xyz = np.array(xyz)

    def __len__(self):
        return len(self.atoms)

    def to_bohr(self):
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.xyz *= 1.889725989

    def to_angstrom(self):
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.xyz /= 1.889725989

    def __str__(self):
        line_form = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
        out = "{:d}\n{:s}\n".format(len(self), self.units)
        for atom, xyz in zip(self.atoms, self.xyz):
            out += line_form.format(atom, *xyz)
        return out 

    def copy(self):
        return Molecule(self.xyz, self.units)

if __name__ == '__main__':
    geom = open('../../extra-files/molecule.xyz','r')
    mol = Molecule(geom)
    mol.to_bohr()
    mol.to_angstrom()
    print(mol)

    
    
