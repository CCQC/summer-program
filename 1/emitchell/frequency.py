#!bin/python

import os
import sys
import numpy as np
import cmath
from numpy import linalg as LA

sys.path.append('../../0/emitchell/')
from molecule import Molecule

class Frequency(object):   # Builds Frequency class
	def __init__(self,hess_file,xyz_file):   # Initializes class with Hessian file specified and .xyz file specified
		with open(hess_file) as f:   # Reads in Hessian file
			hess = [line.split() for line in f]   # Makes Hessian into a list of lists
		try:   # Makes sure Hessian is only numbers
			hess = [list(map(float,i)) for i in hess]   # Makes Hessian into numbers
		except:
			sys.exit("Error: Hessian file must consist of real numbers only.") 	
		self.hess = np.array(hess)
		try:   # Reads in .xyz file and checks if file exists
			self.mol = Molecule(xyz_file)   # Sends .xyz file to molecule.py
		except:
			sys.exit("Error: xyz file given cannot be found.")

	def frequencies(self,Mol,Hess):   # Gets the normal modes and coordinates to be put into a usable JMOL file
		# Constants
		Eh = 4.359784E-18
		BohrR = 5.29177249E-11
		c = 2.99792458E10
		amu_to_kg = 1.6605402E-27 
		# Builds mass matrix
		M = [i for i in Mol.masses for j in range(Mol.natom)]
		M = np.array(M) ** (-1/2)
		M = np.diag(M)
		try:   # Creates the mass-weighted Hessian and checks if the matrices have compatible dimensions
			self.mwHess = LA.multi_dot([M,Hess,M])
		except:
			sys.exit("Error: Dimensions of molecule and Hessian do not match. Check files.")
		vals, vecs = LA.eigh(self.mwHess)   # Solves the mass-weighted Hessian and returns the eigenvalues and eigenvectors
		evecs = M.dot(vecs)   # Builds the un-mass-weighted eigenvectors
		Q = np.transpose(evecs)   # Puts the eigenvectors into columns
		evals = vals * (Eh / ((BohrR ** (2)) * amu_to_kg ))   # Converts the eigenvalues into rad^2 / s^2 units
		k = [cmath.sqrt(i) for i in evals] 
		K = []
		for i in k:   #   Checks if imaginary frequencies are present
			if np.imag(i) != 0:
				i = np.imag(i)*1j   # Makes the imaginary value usable
			else:
				i = np.real(i)
			K.append((1 / (2 * np.pi * c)) * i)   # Converts eigenvalues into wavenumbers
		self.Jmol(K,Q,Mol) 

	def Jmol(self,Evals,Evecs,Mol):   # Writes a usable JMOL file with the values and vectors from frequencies
		f = open('JMOL','w')
		Mol.to_Angstrom()   # Changes coordinates to Angstroms
		string = ""
		for i in range(3 * Mol.natom):    # Loop to create the file
			string += str(Mol.natom) + "\n"
			if np.imag(Evals[i]) != 0:
				string += "{:7.2f}i cm^-1\n".format(np.imag(Evals[i]))
			else:
				string += "{:7.2f} cm^-1\n".format(Evals[i])
			k = 0   # Allows normal coordinates to be written in the same loop as the geometric coordinates
			for j in range(Mol.natom):
				string += "{:<2} {: 14.10f} {: 14.10f} {: 14.10f} {: 14.10f} {: 14.10f} {: 14.10f}\n".format(Mol.labels[j],Mol.geom[j][0],Mol.geom[j][1],Mol.geom[j][2],Evecs[i][k],Evecs[i][k+1],Evecs[i][k+2])
				k += 3
			string += "\n"
		f.write(string)
		f.close()

if __name__ == "__main__":		
	freq = Frequency('../../extra-files/hessian.dat','../../extra-files/molecule.xyz')
	freq.frequencies(freq.mol,freq.hess)
