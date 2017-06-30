# coding: UTF-8

import copy
import numpy
import sys
sys.path.insert(0, '../../extra-files')
import masses

def filename_to_lines(filename):
    '''Takes a filename as input and returns the lines of the file.'''
    try:
        with open(filename, 'r') as file_contents:
            return file_contents.read().splitlines()
    except IOError:
        raise IOError('Error! File {} not found.'.format(filename))

def validate_int(number, err_msg):
    '''Returns the input as an integer if possible and raises an error
    otherwise.'''
    try:
        return int(number)
    except ValueError:
        raise ValueError("{} is not an integer.".format(err_msg))

def get_elt_0(array, target, sought):
    '''Returns the zeroth element of an array if possible and raise an error
    otherwise.'''
    try:
        return array[0]
    except IndexError:
        raise IndexError('{} not found in {}.'.format(target, sought))

def validate_array_length(length, array, name, elt_type):
    '''Verify that the received array is of the proper length and raise an
    error otherwise.'''
    if not isinstance(length, int):
        raise Exception("The array length provided is not an integer.")
    if len(array) != length:
        raise Exception("{} is expected to be {} {}.".format(name, length, elt_type))

def validate_atom(atom, i):
    '''Verify that the received string is an atom symbol and raise an error
    otherwise.'''
    if atom.upper() not in masses.charge:
         raise Exception("Atom {} on line {} not found.".format(atom, i))

def validate_coord(coord, i, j):
    '''Return the received object as a float if possible and raise an error
    otherwise.'''
    try:
        return float(coord)
    except ValueError:
        raise ValueError("Coordinate {} on line {} must be a real number.".format(j, i))

def validate_coord_line(line, i):
    '''Validate the received string as a valid XYZ coordinate line. Return the
    values if possible and raise an error otherwise.'''
    split_line = line.split()
    validate_array_length(4, split_line, "Line {}".format(i), "whitespace separated terms")
    atom_type = split_line.pop(0)
    validate_atom(atom_type, i)
    coords = []
    for j, coord in enumerate(split_line, start=1):
         coords.append(validate_coord(coord, i, j))
    return atom_type, coords

class Molecule(object):
    '''A molecule with atoms of specific masses and charges at specific locations.'''    

    def __init__(self, xyz):
            self.units = "Angstrom"
            self.labels = []
            self.masses = []
            self.charges = []
            xyz_lines = filename_to_lines(xyz)
            # Get the number of atoms.
            natom = get_elt_0(xyz_lines, "Line 1", "file {}".format(xyz)) 
            self.natom = validate_int(natom, "Line 1")
            # Validate the file length.
            validate_array_length(self.natom + 2, xyz_lines, "XYZ input", "lines")
            geom = numpy.empty([0, 3])
            for i, line in enumerate(xyz_lines[2:]):
            	# Ensure a coordinate line is valid before processing it.
                atom, matr_row = validate_coord_line(line, i + 3)
                geom = numpy.vstack((geom, matr_row))
                self.labels.append(atom)
                self.charges.append(masses.get_charge(atom))
                self.masses.append(masses.get_mass(atom))
            self.geom = geom

    def to_bohr(self):
        '''Convert Angstrom units to Bohr units.'''
        if self.units == "Bohr":
            pass
        self.geom = 1.88973 * self.geom
        self.units = "Bohr"
    
    def to_angstrom(self):
        '''Convert Bohr units to Angstrom units.'''
        if self.units == "Angstrom":
            pass
        self.geom = 1 / 1.88973 * self.geom
        self.units = "Angstrom"
    
    def xyz_string(self):
        '''Print the molecular coordinates in xyz format.'''
        for label, geom in zip(self.labels, self.geom.tolist()):
            print("{:3} {:15.10f} {:15.10f} {:15.10f}".format(label, geom[0], geom[1], geom[2]))
    
    def copy(self):
        '''Return a deep copy of the molecule object.'''
        return copy.deepcopy(self)

if __name__ == "__main__":
	water = Molecule("../../extra-files/molecule.xyz")
	water.to_bohr()
	water2 = water.copy()
	water2.to_angstrom()
	water.xyz_string()
	water2.xyz_string()
	print(water.masses)
	print(water.charges)