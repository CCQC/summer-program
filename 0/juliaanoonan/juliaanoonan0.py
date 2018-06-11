import sys
sys.path.insert(0, '../../extra-files')
import masses  #provides functions get_mass and get_charge
import numpy  #lets me make an array later

class Molecule:

    def __init__(self, xyzfile):
        with open(xyzfile, "r") as f:
            lines = f.readlines() #opening the xyz file for reading
        self.xyz = xyzfile

        self.natom = int(lines[0]) #first line = num of atoms
        
        if len(lines) != (self.natom + 2):
            print("Error: incorrectly formatted xyz file")
            return

        self.units = lines[1].rstrip() #second line = Angstroms or Bohr
        if self.units != "Angstrom" and self.units != "Bohr":
            print("Error: unspecified units or units we are currently unequipped to handle.")
            return

        self.labels = []
        for line in lines[2:]:  #from the third line to the end
            self.labels = self.labels + [line[0]] #append the first character
            
        self.masses = [masses.get_mass(i) for i in self.labels] #creates a list of the masses

        self.charges = [masses.get_charge(i) for i in self.labels] #creates a list of the charges

        fgeoms = []
        for line in lines[2:]:  #from the third line to the end
            line = line.split()
            newline = line[1:]  #leave out the atom name, have list of xyz's
            for i in range(len(newline)):
                newline[i] = float(newline[i])  #converting to float
            fgeoms = fgeoms + [newline]  #a list of lists of the xyz values
        self.geoms = numpy.array(fgeoms)  #making an array out of the list of lists

    def to_bohr(self):
        if self.units == "Angstrom":
            self.geoms = self.geoms*(1.8897259886) #conversion to Bohr
            self.units = "Bohr"

    def to_angstrom(self):
        if self.units == "Bohr":
            self.geoms = self.geoms*(0.529177249) #conversion to Angstroms
            self.units = "Angstrom"


    def xyz_string(self):
        with open(self.xyz,"r") as f:
            thestring = f.read()  #reads the file in a string
        return thestring

    def copy(self):
        molecule = Molecule(self.xyz)
        return molecule
        
if __name__ == "__main__":
    water = Molecule("../../extra-files/molecule.xyz")  #create the water
    print(water.units)  #printing the member variables
    print(water.natom)
    print(water.labels)
    print(water.masses)
    print(water.charges)
    print(water.geoms)


    water2 = water.copy()
    water2.to_bohr() #converting to Bohr
    print(water2.units)
    print(water2.geoms)

    water2.to_angstrom() #converting back to Angstroms
    print(water2.units)
    print(water2.geoms)

    print(water.xyz_string()) #printing the xyz string

    water5 = water.copy() #water5 is a copy of water
    water.to_bohr() #changing water's units
    print(water5.units) #but water5 is unchanged!

            
         
