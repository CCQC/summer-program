#__author__ = "mbowman"


import sys
sys.path.insert(0, '../../extra-files')

import numpy

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

if __name__ == "__main__":
	mole_str = open("../../extra-files/molecule.xyz").read()
	mol = Molecule(mole_str)
	print mol.units
