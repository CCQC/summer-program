#! /usr/bin/python3

import numpy as np


class Molecule:
    def __init__(self, geom_str, units='Angstrom'):
        self.read(geom_str)
        self.units = units

    def __len__(self):
        return len(self.atoms)

    def __str__(self):
        line_form = "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n"
        out = "{:d}\n{:s}\n".format(len(self), self.units)
        for atom, coord in zip(self.atoms, self.geom):
            out += line_form.format(atom, coord[0], coord[1], coord[2])
        return out    

    def read(self, geom_str):
        self.atoms = []
        self.geom = []
        with open(geom_str) as coord:
            lines = []
            for line in coord:
                lines.append(line)    
        for info in lines[2:]:
            atom, x, y, z = info.split()
            self.atoms.append(atom)
            self.geom.append([float(x), float(y), float(z)])
        self.geom = np.array(self.geom)

    def angstrom(self):
        if self.units == 'Bohr':
            self.geom = self.geom/1.889725989
            self.units = 'Angstrom'

    def bohr(self):
        if self.units == 'Angstrom':
            self.geom = self.geom*1.889725989
            self.units = 'Bohr'
    


