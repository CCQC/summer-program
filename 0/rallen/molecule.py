
import sys

sys.path.append('../../extra-files')
import numpy as np
from masses import get_charge, get_mass

class Molecule(object):
    def __init__(self, geom_string):
	self.geom_string = geom_string  	
        self.lines = geom_string.strip().split("\n")
        self.natom = int(len(self.lines[2:]))
        self.units = self.lines[1]
        self.labels = []
	self.get_labels()
	self.atoms = []
        self.geom = []
	charge = [get_charge(atom) for atom in self.labels]
	mass = [get_mass(atom) for atom in self.labels]
            
        for line in self.lines[2:]:
            atom, x, y, z = line.split()
            self.atoms.append(atom)
            self.geom.append([float(x), float(y), float(z)])
            
        self.geom = np.array(self.geom)

        # parses file into lines and counts those to find the number of atoms. Then parse
        # that into atoms and coordinates while appending the label and geometry arrays

    def get_labels(self):
	for line in self.lines[2:]:
	    for char in line[0]:
                self.labels.append(char)

	#reads through each line and prints the first character (atom label)
        
    def to_bohr(self):
        if self.units == "Angstrom":
            self.units = "bohr"
            for i in range(0,len(self.geom)):
                for j in range(0,len(self.geom[0])):
                    self.geom[i,j] = self.geom[i,j]*1.88973

        # changes units variable and then scans through geom array through each row and column, changing the value in question
                            
    def to_angstrom(self):
        if self.units == "bohr":
            self.units = "angstrom"
            for i in range(0,len(self.geom)):
                for j in range(0,len(self.geom[0])):
                    self.geom[i,j] = self.geom[i,j]*0.529177
        # changes units variable and then scans geom array through each row and column, changing the value in question

    def xyz_string(self):
	var = ""
	for i in range(self.natom):
	    var += str(self.labels[i]) + ' ' + '\t'.join([str(j) for j in self.geom[i]]) + "\n"

	return str(self.natom) + '\n' + str(self.units)  + '\n' + var
        # returns the number of atoms in the molecule, then the units, then the label, mass, charge, and x,y,z coordinates of each atom in the molecule

    def copy(self):
        return Molecule(self.xyz_string())

# returns geometry and units from input file

string = open("../../extra-files/molecule.xyz", "r").read()
molecule1 = Molecule(string)
print molecule1.copy().xyz_string()

# tests codE



