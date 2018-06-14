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
		geom = []
		self.labels = []
		self.masses = []
		self.charges = []
		lines = geom_str.strip().splitlines()
		self.natom = int(lines[0])
		for line in lines[2:]:
			atom, x, y, z = line.split() 
			self.labels.append(atom) 
			geom.append([float(x), float(y), float(z)])
		self.geom = np.array(geom)
		for x in self.labels:
			self.masses.append(get_mass(x))
			self.charges.append(get_charge(x))
		y = 0
		for x in self.labels:
			print("Mass of " + x + " is: " + str(self.masses[y]))
			print("Charge of " + x + " is: " + str(self.charges[y]))
			y += 1
#get_mass and get_charge functions not being used (should append into self.masses and self.charges

					

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
		
	
if __name__ == "__main__":	
	file1 = open("../../extra-files/molecule.xyz").read()
	mol = Molecule(file1, "Angstrom")
	print("Units are currently in: " + mol.units)
	print(mol.geom)
	print("Number of atoms: " + str(mol.natom))
	print("Atoms: " + str(mol.labels))
	mol.to_bohr()
	print("Units have now been converted to: " + mol.units)
	print(mol.geom)
	mol.to_angstrom()
	print("Units have now been converted back to: " + mol.units)
	print(mol.geom)
