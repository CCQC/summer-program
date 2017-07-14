import numpy as np
import masses

class Molecule:
    def __init__(self,filename):
        self.moleculefile = open(filename,'r')
        self.moleculelines = self.moleculefile.readlines()
        self.units = str(self.moleculelines[1]).rstrip()

        self.natom = int((self.moleculelines[0]).rstrip())
        self.labels = [0]*self.natom
        self.masses = [int(0)]*self.natom
        self.charges = [float(0)]*self.natom
        self.geom = np.zeros(shape=(self.natom,3))
        self.geom[1][0] = 5

        for atom in range(0,self.natom):
            tempatom = self.moleculelines[atom+2].split()
            print(tempatom)
            self.labels[atom] = str(tempatom[0])
            self.charges[atom] = int(masses.get_charge(self.labels[atom]))
            self.masses[atom] = float(masses.get_mass(self.labels[atom]))

            for i in range(1, 4):
                self.geom[atom][i-1] = tempatom[i]
    def to_bohr(self):
        if self.units == "Bohr":
            return 0
        
        elif self.units == "Angstrom":
            self.units = "Bohr"
            self.geom *= 1.889725989
            return 0
    def to_angstrom(self):
        if self.units == "Angstrom":
            return 0
        elif self.units == "Bohr":
            self.units = "Angstrom"
            self.geom *= (1./1.889725989)
    def xyz_string(self):
        outstring = str(self.natom) + "\n" + str(self.units) + "\n" + \
            str(self.geom)
    def copy(self):
        return self
