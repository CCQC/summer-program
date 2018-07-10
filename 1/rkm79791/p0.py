#This is the first code I've ever written in python... maybe it works! 
#By Ryan Maynard
import copy
import numpy
from masses import get_charge, get_mass
class Molecule:
    def __init__(self, file_name):
        self.file_name = file_name
        with open(file_name) as geom_file:
            info = geom_file.readlines()
        info = [x.strip() for x in info] #I wanted to strip the \n but keep it a list
        #print(info) 
        self.info = info


        self.units = info[1]
        
        self.num_atoms = info[0]
        
        #empty lists that I append to in the code below them
        atom_labels = []        
        geom = []
        #here I am breaking up the "pieces" inside the list "info" which is a list consisting of each line of the xyz file
        for piece in info[2:]:
            if not piece == " ":
                piece = piece.split()
                #This appends the first string (the actual atom label) to the list "atom_label"
                atom_labels.append(piece[0])
                #I made a "temporary list" to help myself think through the process of converting the string coordinates in the mini-list to floats so that they can be used later in out to_bohr and to_angstrom functions 
                temp_list = [float(x) for x in piece[1:]]
                #I then appended to the fittingly-named list, "geom"
                geom.append(temp_list)
            
        self.atom_labels = atom_labels
        self.geom = numpy.array(geom)
        
        masses = list(map(get_mass, atom_labels))
        self.mass_list = masses

        charges = list(map(get_charge, atom_labels))
        self.charge_list = charges

    #reads the units of the file and converts from angstrom to bohr if necessary
    def to_bohr(self):
        if self.units == "angstrom" or "Angstrom" or "Angstroms" or "angstroms":
            self.geom = self.geom*1.8897259886
            self.units = "Bohr"


    #opposite of to_bohr
    def to_angstrom(self):
        if self.geom == "bohr" or "Bohr" or "Bohrs" or "bohrs":
            angstrom_geom = self.geom/1.8897259886
            self.geom = angstrom_geom
            self.units = "Angstrom"
        else:
            pass
        

#this just returns the initial string with "new line" stripped
    def xyz_string(self):
        New_String = str(self.num_atoms) + "\n" + str(self.units) + "\n" + str(self.geom)
        return New_String

    def copy(self): 
        new_copy = copy.deepcopy(self)
        return new_copy



if __name__ == "__main__":

#I added some spaces and labels because it made me able to think through the process better and ensure the correct data was printing where desired 
    my_mol = Molecule("molecule.xyz")
    print("Units = " + my_mol.units)
    print( )
    print("Number of Atoms = " + my_mol.num_atoms)
    print()
    print("Atoms in Molecule:")
    print(my_mol.atom_labels)
    print( )
    print("List of Masses in Molecule:")
    print(my_mol.mass_list)
    print( )
    print("List of Charges in Molecule:")
    print(my_mol.charge_list)
    print(                               )
    print("Geometry of Molecule (X, Y, Z):")
    print(my_mol.geom)
    my_mol.xyz_string()
#    my_mol.copy()
