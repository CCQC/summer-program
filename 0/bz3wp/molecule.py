# Defines a class molecule which takes an input file and stores the various components.

import sys
sys.path.insert(0,'../../extra-files')
import masses
import numpy as np

class Molecule(object):
    
    def __init__(self, geom_str):
        lines=geom_str.splitlines()                             # splitting lines of input file
        self.natom = int(lines[0])                              # number of atoms in the input
        self.units = lines[1]                                   # units of geometry
        self.labels = []                                        # atoms in the input file
        self.mass = []                                          # mass of each atom
        self.charge = []                                        # charge of each atom
        self.geom = []                                          # geometry of molecule in xyz coordinates
        for line in lines[2:]:
            atom,x,y,z = line.split()
            self.labels.append(str(atom))
            self.mass.append(float(masses.get_mass(atom)))
            self.charge.append(int(masses.get_charge(atom)))
            self.geom.append([float(x),float(y),float(z)])
        self.geom = np.array(self.geom)                         # matrix of xyz coordinates
    
    def to_bohr(self):                                          # function to convert units to Bohr
        if self.units == "Bohr":
            print("Already in Bohr")
        else:
            self.geom *= 1.8897 
            self.units = "Angstrom"
  
    def to_angstrom(self):                                      # function to convert units to Angstrom
        if self.units == "Angstrom":
            print("Already in Angstrom")
        else:
            self.geom *= 0.5292
            self.units = "Bohr"

    def __str__(self):                                          # returns a .xyz file with contents of input file
        out = "{}\n\n".format(self.natom)
        for atom,line in zip(self.labels, self.geom):
            x,y,z=line
            out += "{:2s}{:>10.5f}{:>10.5f}{:>10.5f}\n".format(atom,x,y,z)
        return out
    
    def copy(self):                                             # returns molecule with a copy of itself
        return Molecule(self)
    
    @staticmethod                                               # function to read in input
    def from_file(infile = "geom.xyz"):
        with open(infile) as f:
            return Molecule(f.read())

       
