import numpy as np
import masses

class Molecule:
    """
    defines the Molecule class. Molecule objects are initialized with a .xyz
    filename as input. 
    """
    def __init__(self,filename):
        self.moleculefile = open(filename,'r')
        self.moleculelines = self.moleculefile.readlines()
        
        #units are supplied on line 2 
	self.units = str(self.moleculelines[1]).rstrip()
        #number of atoms is supplied on line 1   
        self.natom = int((self.moleculelines[0]).rstrip())
        
        #create empty lists of atom labels, masses, charge, and geometry matrix         for later population
        self.labels = [0]*self.natom
        self.masses = [int(0)]*self.natom
        self.charges = [float(0)]*self.natom
        self.geom = np.zeros(shape=(self.natom,3))
        
        #loop through atoms and populate information into empty lists
        for atom in range(0,self.natom):
            #+2 because atoms start at line 3
            tempatom = self.moleculelines[atom+2].split()
            
            self.labels[atom] = str(tempatom[0])
            self.charges[atom] = int(masses.get_charge(self.labels[atom]))
            self.masses[atom] = float(masses.get_mass(self.labels[atom]))

            for i in range(1, 4):
                #-1 for indexing reasons
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
