#Author: Ally Burke


import sys
sys.path.append('../../extra-files')
import copy
import numpy as np
from masses import get_charge
from masses import get_mass

class Molecule(object):

	#Used to initialize an instance of the Molecule class

	def __init__(self, geom_str, units):
		self.units = units				#Specifies units for spatial coordinates
		self.read(geom_str)				#Reads in the geometry file for the specific molecule		

	#Returns a string representing the contents of the Molecule object in .xyz format
	
	def read(self, geom_str):
		#Read in the file containing the geometry of the molecule in question
		print("\n" + "MOLECULE INFORMATION" + "\n")
		geom = []
		self.labels = []
		self.masses = []
		self.charges = []
		lines = geom_str.strip().splitlines()
		self.natom = int(lines[0])
		print("Number of atoms: " + str(self.natom) + "\n")
		print("Units are currently in: " + self.units + "\n")
		for line in lines[2:]:
			atom, x, y, z = line.split() 
			self.labels.append(atom) 
			geom.append([float(x), float(y), float(z)])
		self.geom = np.array(geom)	
		print("Atoms: " + str(self.labels)+ "\n")
		
		#Reading in the masses and charges of the atoms from the .xyz file
		
		for x in self.labels:
			self.masses.append(get_mass(x))
			self.charges.append(get_charge(x))
		y = 0
		for x in self.labels:
			print("Mass of " + x + " is: " + str(self.masses[y]))
			print("Charge of " + x + " is: " + str(self.charges[y]))
			y += 1
		print("\n")
		print("Molecular Geometry in: " + self.units + "\n")
		print(self.geom)


	#UNIT CONVERSION FUNCTIONS

	def to_bohr(self):

	#Converts units from Angstrom to Bohr	
		
		if self.units == "Angstrom":
			self.units = "Bohr"
			self.geom *= 1.889725989
			
	def to_angstrom(self):
		
	#Converts from Bohr to Angstrom

		if self.units == "Bohr":
			self.units = "Angstrom"
			self.geom /= 1.889725989

	#Create fresh copy of Molecule object
	
	def copy(self):
		x = Molecule(self.read(geom_str))
		return x			
		
#Used to test the functions above

if __name__ == "__main__":	
	file1 = open("../../extra-files/molecule.xyz").read()
	mol = Molecule(file1, "Angstrom")
	
	#Testing unit conversion functions
	mol.to_bohr()
	print("\n" + "Units have now been converted to: " + mol.units)
	print(mol.geom)
	mol.to_angstrom()
	print("\n" + "Units have now been converted back to: " + mol.units)
	print(mol.geom)
