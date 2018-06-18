#!/usr/bin/env python
import sys 
import numpy as np

sys.path.insert(0, '../../extra-files')
from masses import get_charge, get_mass


class Molecule(object):
    
    def __init__(self,geom_str,units = "Angstrom"):
    
        self.read(geom_str)
        self.units = units
        self.charge = [get_charge(i) for i in self.label]
        self.mass   = [get_mass(i)   for i in self.label]

    def read(self,geom_str):
        geom = []
        self.label = []
        lines = geom_str.strip().splitlines()
        natom = int(lines[0])
        for line in lines[2:]:
            atom, x, y, z = line.split() 
            self.label.append(atom.upper())
            geom.append([float(x), float(y), float(z)])
        self.geom = np.array(geom)

    def to_angstrom(self):
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.geom *= 0.5291772085936

    def to_bohr(self):
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.geom /= 0.5291772085936

    def __len__(self):
        return len(self.geom)

    def copy(self):
        return Molecule(self.geom, self.units)

    def __str__(self):
        line_form = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
        out = "{:d}\n{:s}\n".format(len(self), self.units)
        for atom, xyz in zip(self.label, self.geom):
            out += line_form.format(atom, *xyz)
        return out

if __name__ == '__main__':
    geom_str = open("../../extra-files/molecule.xyz","r").read()
    mol = Molecule(geom_str,"Angstrom")
    mol.to_bohr()
print(mol)


