#__author__ = "mbowman"


import sys
sys.path.insert(0, '../../extra-files')

import masses as m

import numpy as np

class Molecule(object):
 
	def __init__(self, mole_str, units="angstrom"):
		"""
		Generates an instance of a molecule
		:param 	
		"""	
		self.parse(mole_str)

	def parse(self, mole_str):
		lines = mole_str.strip().split("\n")
		self.natom = lines[0]
		if lines[1].lower() == "angstrom":
			self.units = "angstrom"
		elif lines[1].lower() == "bohr":
			self.units = "bohr" 
		self.labels = [i.strip().split(" ")[0].strip() for i in lines[2:]]
		self.masses =  [m.get_mass( l ) for l in self.labels ]
		self.charges = [m.get_charge( l ) for l in self.labels ]
		tempArray = [i.split("   ")[1:4] for i in lines[2:]]
		for a in range(len(tempArray)):
			for b in range(len(tempArray[a])):
				tempArray[a][b] = float(tempArray[a][b]
)
		self.geom = np.array(tempArray)

if __name__ == "__main__":
	mole_str = open("../../extra-files/molecule.xyz").read()
	mol = Molecule(mole_str)
	print mol.units
	print mol.labels
	print mol.masses
	print mol.charges
	print mol.geom
