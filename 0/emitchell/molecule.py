#!bin/python

import os
import sys
import numpy as np

sys.path.append('../../extra-files/')
import masses

class Molecule(object):
	def __init__(self,xyz_file): #Initializes member variables: natom, units, labels, masses, charges, geom
		with open(xyz_file) as f:
			geom = [line.split() for line in f]
		try: #Tests formatting of .xyz file
			self.natom = int(geom.pop(0)[0]) #Gets number of atoms
			units = geom.pop(0) #Gets units of coordinates
			for i in units:
				if isinstance(i, str) == True:
					i = str(i)
					if i.upper() == "BOHR" or i.upper() == "ANGSTROM":
						self.units = i
		except:
			sys.exit("Error: Check .xyz file to ensure correct format. \n Number of Atoms \n Unit of Length \n Atom with Coordinates \n . \n . \n .")
		self.labels = []
		self.masses = []
		self.charges = []
		for i in geom: #Loops over the labels with coordinates of the file
			a = i.pop(0) 
			self.labels.append(a) #Builds list of labels
			try:
				self.masses.append(masses.get_mass(a)) #Builds list of atomic masses
				self.charges.append(masses.get_charge(a)) #Builds list of atom charges
			except:
				sys.exit("Error: Check that labels are written correctly.")
		try:
			geom = [list(map(float,j)) for j in geom] 
			self.geom = np.array(geom) #Makes array out of list of lists containing coordinates
		except:
			sys.exit("Error: Check if coordinates are in numerical format with x, y, and z components.")
		f.close()
	
	def to_Bohr(self): #Converts Angstrom to Bohr
		if self.units.upper() == "ANGSTROM": #Checks if in Angstroms
			self.geom *= 1.88973
			self.units = "Bohr"
		elif self.units.upper() != "BOHR": #Checks if in usable units
			print("Error: Units must be in Angstrom or Bohr")
		else:
			print("Already in Bohr")
		return self

	def to_Angstrom(self): #Converts Bohr to Angstrom
		if self.units.upper() == "BOHR": #Checks if in Bohr
			self.geom *= 0.529177 
			self.units = "Angstrom"
		elif self.units.upper() != "ANGSTROM": #Checks if in usable units
			print("Units must be in Angstrom or Bohr")
		else:
			print("Already in Angstrom")
		return self
	
	def xyz_string(self): #Gives a .xyz file
		string = "units " +  str(self.units.lower())
		for i in range(self.natom):
			string += "\n {:^2} {:^20.10f} {:^20.10f} {:^20.10f}".format(self.labels[i],self.geom[i][0],self.geom[i][1],self.geom[i][2])
		return string 
	
	def copy(self): #Returns fresh copy of self
		f = open('copy.xyz','w')
		f.write(repr(self))
		f.close()
		return Molecule('copy.xyz')
	
	def __len__(self): #Allows the Molecule to be used by len()
		return self.natom
	
	def __str__(self): #Allows the Molecule to be turned into a string
		s = "Number of atoms: " + str(self.natom) + "\n" + "Unit of Length: " + str(self.units) + "\n"
		s += "\n {:^10} | {:^15} | {:^10} | {:^20} | {:^20} | {:^20}".format("Label","Mass","Charge","X-Coordinate","Y-Coordinate","Z-Coordinate")
		s += "\n" + ("-" * 110)
		for i in range(self.natom):
			s += "\n {:^10} | {:^ 15.10f} | {:^ 10.2f} | {:^ 20.10f} | {:^ 20.10f} | {:^ 20.10f}".format(self.labels[i],self.masses[i],self.charges[i],self.geom[i][0],self.geom[i][1],self.geom[i][2])
		return s

	def __repr__(self): #Allows the Molecule to be used by repr()
		string = str(self.natom) + "\n" + str(self.units)
		for i in range(self.natom):
			string += "\n {:^2} {:^20.10f} {:^20.10f} {:^20.10f}".format(self.labels[i],self.geom[i][0],self.geom[i][1],self.geom[i][2])
		return string 

	def __iter__(self): #Initializes the for loop
		self.n = 0
		myList = []
		for i in range(self.natom):
			item = (self.labels[i],self.geom[i])
			myList.append(item)
		self.myTuple = tuple(myList)
		return self

	def __next__(self): #Iterates through the conditions set by __iter__ when used in for loop
		if self.n < self.natom:
			line = self.myTuple[self.n]
			self.n += 1
			return line
		else:
			raise StopIteration

	def __add__(self,other): #Let's two .xyz files be combined
		self.natom += other.natom
		if self.units != other.units:
			sys.exit("Error: Files must have same units to be combined.")
		self.labels = self.labels + other.labels
		for i, j in enumerate(other.geom):
                        for l, k in enumerate(j):
                                other.geom[i][l] = str(float(k) + 1)
		self.geom = np.concatenate((self.geom,other.geom), axis=0)
		self.masses += other.masses
		self.charges += other.charges
		return self
	

if __name__ == "__main__":
	mol = Molecule('../../extra-files/molecule.xyz')
	print("Converting from Angstrom to Bohr: ")
	mol.to_Bohr()
	print(repr(mol) + "\n")
#	print("Converting back to Angstrom: ")
#	mol.to_Angstrom()
	print(repr(mol) + "\n")
	print("Copy of self saved in copy.xyz: ")
	mol.copy()
	print(repr(mol) + "\n")
	print("How many atoms are in our molecule?: " + str(len(mol)) + "\n")
	print("Pretty way to read contents of our molecule: ")
	print(str(mol) + "\n")
	print("Let's iterate through: ")
	for i in mol:
		print(i)
	print("\n")
	print("And lastly, we can add two molecules: ")
	mol2 = Molecule('copy.xyz')
	mol + mol2
	print(mol)	
