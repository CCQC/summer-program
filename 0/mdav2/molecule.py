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
        #checking if input string is filename or xyzstring by determining number of lines
        if data.endswith('.xyz'):
            with open(data, 'r') as moleculefile:
                self.moleculelines = moleculefile.readlines()
        else:
            self.moleculelines = data.split('\n')
        #units are supplied on line 2
        self.units = self.moleculelines[1].rstrip()
        #number of atoms is supplied on line 1
        self.natom = int((self.moleculelines[0]).rstrip())

        #create empty lists of atom labels, masses, charge, and geometry matrix         for later population
        self.labels = []
        self.masses = []
        self.charges = []
        self.geom = np.zeros((self.natom,3))

        #loop through atoms and populate information into empty lists
        for atom in range(self.natom):
            #+2 because atoms start at line 3
            element,x,y,z = self.moleculelines[atom+2].split()

            self.labels.append(element)
            self.charges.append(masses.get_charge(element))
            self.masses.append(masses.get_mass(element))
            self.geom[atom] = [x,y,z]

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
        #returns a string containing the object in .xyz format.
        outstring = '{}\n{}\n'.format(self.natom,self.units)
        
        for atom in range(0, self.natom):
            outlist = [self.labels[atom]] + self.geom[atom].tolist()
            outstring += ('{:6}' + '{: 012.11f}   '*3+'\n').format(self.labels[atom],*self.geom[atom])
        return outstring

    def copy(self):
        return Molecule(self.xyz_string())
