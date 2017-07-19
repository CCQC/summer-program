import sys
sys.path.append('../../extra-files')
import numpy as np
import masses

class Molecule:
    """
    defines the Molecule class. Molecule objects are initialized with a .xyz
    filename as input.
    """
    def __init__(self,data):
        #checking if input string is filename or xyzstring
        if len(data.split('\n')) != 1:
            self.moleculelines = data.split('\n')
        else:
            with open(data, 'r') as moleculefile:
                self.moleculelines = moleculefile.readlines()

        #units are supplied on line 2
        self.units = str(self.moleculelines[1]).rstrip()
        #number of atoms is supplied on line 1
        self.natom = int((self.moleculelines[0]).rstrip())

        #create empty lists of atom labels, masses, charge, and geometry matrix         for later population
        self.labels = [0]*self.natom
        self.masses = [0]*self.natom
        self.charges = [0.0]*self.natom
        self.geom = np.zeros((self.natom,3))

        #loop through atoms and populate information into empty lists
        for atom in range(0,self.natom):
            #+2 because atoms start at line 3
            element,x,y,z = self.moleculelines[atom+2].split()

            self.labels[atom] = str(element)
            self.charges[atom] = int(masses.get_charge(self.labels[atom]))
            self.masses[atom] = float(masses.get_mass(self.labels[atom]))

            for i in range(0, 3):
                #-1 for indexing reasons
                self.geom[atom][i] = [x,y,z][i]
    def to_bohr(self):
        #converts all units to Bohr radii when called
        #1 br = 0.52917724888
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.geom *= 1.889725989

    def to_angstrom(self):
        #converts all units to Angstrom's when called
        #1 A = 1.889725989 br
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.geom *= (1./1.889725989)
    def xyz_string(self):
        #returns a string
        outstring = ""
        outstring += str(self.natom) + "\n"
        outstring += str(self.units) + "\n"
        for atom in range(0, self.natom):
            outstring += str(self.labels[atom]) + "     "
            outstring += '{: 012.11f}'.format(self.geom[atom][0]) + "   "
            outstring += '{: 012.11f}'.format(self.geom[atom][1]) + "   "
            outstring += '{: 012.11f}'.format(self.geom[atom][2]) + "   "
            outstring += "\n"
        return outstring
    def copy(self):
        return Molecule(self.xyz_string())
