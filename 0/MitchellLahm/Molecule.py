import masses
import numpy as np
from copy import deepcopy

"""This class will store the geometry, masses, and nuclear charges of a molecule."""
class Molecule(object):
	def __init__(self,TextFile):
		# Empty arrays defined to be later appended to
		self.textfile = TextFile
		labels = []
		geom = []
		charges = []
		mass = []
		# open text file, read in the values, and store them as the desired units
		with open(TextFile,'r') as my_file:
			data = my_file.read().replace('\n',' ')
			data = data.split()
			self.natom = int(data[0])
			del data[0]
			self.units = data[0]
			del data[0]
			# reading in rest of the data
			for i in range(self.natom):
				# Store the first line of the column as a label
				labels.append(data[4*i])
				# Temporary list defined to store subarray of geom
				temp = []
				# Cycle through the x,y,z coordinates
				for j in range(3):
					temp.append(float(data[4 * i + j + 1]))
				# Append 1x3 subarray to geom array
				geom.append(temp)
			self.labels = labels
			self.geom = np.array(geom)
			for i in range(self.natom):
				mass.append(masses.get_mass(str(labels[i])))
				charges.append(masses.get_charge(str(labels[i])))
			self.charges = charges
			self.masses = mass
	
	# Method to convert the units of the cartesian coordinates from Angstrom to Bohr
	def to_bohr(self):
		if self.units != 'Bohr':
			self.units = 'Bohr'
			for i in range(self.natom):
				for j in range(self.natom):
					self.geom[i][j] *= 1.88971616463
	
	# Method to convert the units of the cartesian coordinates from Bohr to Angstrom
	def to_angstrom(self):
		if self.units != 'Angstrom':
			self.units = 'Angstrom'
			for i in range(self.natom):
				for j in range(self.natom):
					self.geom[i][j] /= 1.88971616463
	
	# Method to print the acquired values in the same format as the .xyz file
	def xyz_string(self):
		print str(self.natom)
		print str(self.units)
		for i in range(self.natom):
			x = self.geom[i]
			print '%s    %s' % (self.labels[i],'   '.join('%13.10f' % j for j in x))
	
	# Method to make a fresh copy of Molecule whose values can be changed independently
	# of the initial Molecule class.		
	def copy(self):
		molecule = Molecule(self.textfile)
		molecule.units = deepcopy(self.units)
		molecule.natom = deepcopy(self.natom)
		molecule.labels = deepcopy(self.labels)
		molecule.masses = deepcopy(self.masses)
		molecule.charges = deepcopy(self.charges)
		molecule.geom = deepcopy(self.geom)
		return molecule

# Defining first molecule 'm'
m = Molecule('molecule.xyz')
print m.geom
print m.units

# Converting the units of m from Angstrom to Bohr
m.to_bohr()
print m.geom
print m.units

# Fresh copy of m made, stored as 'n'
n = m.copy()
print n.geom
print n.units

# m converted back to Angstrom, but n retains Bohr units
m.to_angstrom()
print m.geom
print m.units
print n.geom
print n.units

# xyz string printed in same format as the supplied xyz file.
m.xyz_string()
