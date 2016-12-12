#!/usr/bin/env python

import re,regex,sys,numpy as np
from masses import get_mass, get_charge

class Molecule:
	"""
	Object that contains the basic properties of a molecule
	Employs internal coordinate representation
	"""

	def __init__(self,inputString,lengthUnits="Angstrom",angleUnits="Degree"):

		self.inputString = inputString
		self.lengthUnits = lengthUnits
		self.angleUnits  = angleUnits

		self.read(inputString)

		self.masses  = [ float(get_mass(i)) for i in self.atoms ]
		self.charges = [ int(get_charge(i)) for i in self.atoms ]


	def read(self,inputString):

		inputLines     = [i.strip() for i in inputString.splitlines()[1:] if i.strip() != ""]

		self.atoms      = []

		for line in inputLines:
			firstEntry = line.split()[0]
			if regex.isAtom(firstEntry):
				self.atoms.append(firstEntry)
		
		self.Natom  = len(self.atoms)

		inputString  = inputString.replace('{', '{{').replace('}', '}}')

		match = None
		for match in re.finditer(regex.Zmatrix(), inputString, re.MULTILINE):
			pass
		if not match:
			raise Exception('Cannot find start of internal coordinate block')

		start, end       = match.start(), match.end()
		zmatBlock        = inputString[start:end]
		header,footer    = inputString[:start],inputString[end:]
		bodyTemplate     = re.sub('\s+(-?\d+\.?\d+)', ' {:> 17.12f}', zmatBlock)

		coords           = [float(i) for i in re.findall(regex.coordinate,zmatBlock)]
		self.coords      = np.array(coords)
		self.labels      = [i for i in re.findall(regex.coordLabel,zmatBlock)]
		self.Template    = header + bodyTemplate + footer

		
	def __len__(self):
		return self.Natom


	def __str__(self):
		return self.Template.format(*self.coords)
		

	def toBohr(self):
		"""
		:return: 1darray of internal coordinates with bond distances in bohr and angles in radians
		"""
		N = self.Natom

		if self.lengthUnits == "Angstrom":
			for i in range(N-1):
				self.coords[i] /= 0.529177
			self.lengthUnits = "Bohr"

		return self.coords

	def toAngstrom(self):
		
		N = self.Natom
		if self.lengthUnits == "Bohr":
			for i in range(N-1):
				self.coords[i] *= 0.529177
			self.lengthUnits = "Angstrom"

		return self.coords

	def toRadian(self):
		
		N = self.Natom
		if self.angleUnits == "Degree":
			for i in range(N-1,3*N-6):
				self.coords[i] *= np.pi/180
			self.angleUnits = "Radian"

		return self.coords

	def toDegree(self):

		N = self.Natom
		if self.angleUnits == "Radian":
			for i in range(N-1,3*N-6):
				self.coords[i] *= 180/np.pi
			self.angleUnits = "Degree"

		return self.coords

	def copy(self):
		"""
		:return: a fresh copy of the molecule object
		"""
		return Molecule(str(self),self.lengthUnits,self.angleUnits)


if __name__ == '__main__':
	f   = open("h2o.zmat","r").read()
	mol = Molecule(f)
	mol.toBohr()
	print( str(mol) )
