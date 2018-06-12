#Author: Ally Burke


import sys
sys.path.append('../../extra-files')
import copy
import numpy as np
from masses import get_charge
from masses import get_mass

class Molecule(object):

	#Used to initialize an instance of the Molecule class

	def __init__(self, geom_str, units = "Angstroms"):
		self.units = units				#Specifies units for spatial coordinates
		self.read(geom_str)				#Reads in the geometry file for the specific molecule	
		self.natom = self.read(geom_str)		#Specifies the number of atoms in the molecule
		self.labels = self.read(geom_str)		#Specifies list of atom symbols following the input file	
		
		#unsure how to reference parameters without hardcoding
		self.masses = get_mass("H") 	        	#Reads in masses from input file
		self.charges = get_charge("H")			#Reads in charges from input file

	#Unsure if read function is necessary for this project

	def read(self, geom_str):
		
		#Read in the file containing the geometry of the molecule in question
	
        	geom = []
        	self.label = []
        	lines = geom_str.strip().splitlines()
        	natom = int(lines[0])
        	for line in lines[2:]:
           		atom, x, y, z = line.split() 
           		self.label.append(atom)
           		geom.append([float(x), float(y), float(z)])
        	self.geom = np.array(geom)

	#Returns a string representing the contents of the Molecule object in .xyz format

	def xyz_string(self):
		string = "{}\n{}\n".format(self.natom, self.units)
		lines = geom_str.strip().splitlines()
		for line in lines:
			atom, x, y, z = line.split()
			string.append(atom)
			string.append("   " + x + "   " + y + "   " + z)
		return string
					

	#UNIT CONVERSION FUNCTIONS

	def to_bohr(self):

	#Converts units from Angstrom to Bohr	

		if self.units == "Angstrom":
			self.units = "Bohr"
			self.units *= 1.889725989
		elif self.units == "Bohr":
			print("Units already in Bohr")

	def to_angstrom(self):
		
	#Converts from Bohr to Angstrom

		if self.units == "Bohr":
			self.units = "Angstrom"
			self.units /= 1.889725989
		elif self.units == "Angstrom":
			print("Units already in Angstrom")

	#Create fresh copy of Molecule object
	
	def copy(self):
		x = new_geom_str.copy(geom_str)
		return x			
		
	
if __name__ == "__main__":	
	file1 = open("../../extra-files/molecule.xyz").read()
	mol = Molecule(file1)
	print(mol)
