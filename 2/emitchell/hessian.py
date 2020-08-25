#!bin/python

import os
import sys
import re
import time
import subprocess
import numpy as np

sys.path.append('../../0/emitchell/')
from molecule import Molecule
sys.path.append('../../extra-files/')
sys.path.append('../../1/emitchell/')
from frequency import Frequency

class Hessian(object):
	def __init__(self,xyz_file):   # Initializes molecule and locks in original coordinates
		self.mol = Molecule(xyz_file)	
		self.fyle = xyz_file

	def generate_inputs(self,temp,disp_size,directory):   # Makes a directory of directories with input files for the finite differences
		if self.mol.units.upper() == 'ANGSTROM': 
			disp_size *= 0.529177   # Changes displacement size to be in angstroms
		# No displacements
		path = os.path.join(directory + '/XXX/')
		os.makedirs(path)
		with open('{:s}input.dat'.format(path),'w') as f:
			f.write(temp.format(self.mol.xyz_string()))
		# Diagonal, single displacements
		for i in range(3 * self.mol.natom):
			# Positive displacements
			folder = str(i) + '+/'
			self.create_coords(i,None,temp,directory,folder,disp_size,True,False,False)
			# Negative displacements
			folder = str(i) + '-/'
			self.create_coords(i,None,temp,directory,folder,disp_size,False,True,False)
		# Off-Diagonal, double displacements
		for i in range(3 * self.mol.natom):
			for j in range(i+1,3 * self.mol.natom):
				# Positive - positive displacements
				folder = str(i) + str(j) + '++/'
				self.create_coords(i,j,temp,directory,folder,disp_size,True,False,True)
				# Negative - negative displacements
				folder = str(i) + str(j) + '--/'
				self.create_coords(i,j,temp,directory,folder,disp_size,False,True,True)
	
	def create_coords(self,index1,index2,temp,directory,folder,disp,positive,negative,double):   # Creates the displaced coordinates
		coords = Molecule(self.fyle)   # Calls a fresh copy of the coordinates every iteration
		coords.geom = coords.geom.flatten()
		path = os.path.join(directory + '/' + folder)
		os.makedirs(path)
		if positive == True and double == False:   # Positive single displacement
			coords.geom[index1] += disp
		if negative == True and double == False:   # Negative single displacement
			coords.geom[index1] -= disp
		if positive == True and double == True:   # Positive double displacements
			coords.geom[index1] += disp
			coords.geom[index2] += disp
		if negative == True and double == True:   # Negative double displacements
			coords.geom[index1] -= disp
			coords.geom[index2] -= disp
		coords.geom = np.reshape(coords.geom,(3,3))
		with open('{:s}input.dat'.format(path),'w') as f:
			f.write(temp.format(coords.xyz_string()))
	
	def run_jobs(self,command,directory,sec):   # Goes into each directory to submit the input files
		cmd = command + ' input.dat'
		# No displacements 
		os.chdir(directory + '/XXX/')
		subprocess.call(cmd,shell=True)	
		# Single displacements
		for i in range(3 * self.mol.natom):
			path = '/' + str(i) + '+/'
			os.chdir('../../' + directory + path)
			subprocess.call(cmd,shell=True)	
			path = '/' + str(i) + '-/'
			os.chdir('../../' + directory + path)
			subprocess.call(cmd,shell=True)	
		# Double displacements
		for i in range(3 * self.mol.natom):
			for j in range(i+1,3 * self.mol.natom):
				path = '/' + str(i) + str(j) + '++/'
				os.chdir('../../' + directory + path)
				subprocess.call(cmd,shell=True)	
				path = '/' + str(i) + str(j) + '--/'
				os.chdir('../../' + directory + path)
				subprocess.call(cmd,shell=True)	
		os.chdir('../../')
		time.sleep(sec)   # Sets a buffer time if submitted to a queue

	def build_hessian(self,energy,h,directory):   # Creates the Hessian from the energies of the displaced geometries
		size = 3 * self.mol.natom
		hess = np.asmatrix(np.zeros((size,size)))
		const = ( 1 / (2 * (h ** 2))) 
		# No displacements
		path = directory + '/XXX/output.dat'
		exxx = self.awk(path,energy) 
		# Off diagonal
		for i in range(size):
			for j in range(i+1,size):
				path = directory + '/' + str(i) + str(j) + '++/output.dat'				
				epp = self.awk(path,energy)
				path = directory + '/' + str(i) + str(j) + '--/output.dat'				
				enn = self.awk(path,energy)
				path = directory + '/' + str(i) + '+/output.dat'
				epx = self.awk(path,energy)
				path = directory + '/' + str(i) + '-/output.dat'
				enx = self.awk(path,energy)
				path = directory + '/' + str(j) + '+/output.dat'
				exp = self.awk(path,energy)
				path = directory + '/' + str(j) + '-/output.dat'
				exn = self.awk(path,energy)
				hess[i,j] = const * (epp + enn - epx - enx - exp - exn + (2 * exxx))
		hess = hess + hess.transpose()   # Hessians are symmetric
		# Diagonal
		for i in range(size):
			path = directory + '/' + str(i) + '+/output.dat'
			ep = self.awk(path,energy)
			path = directory + '/' + str(i) + '-/output.dat'
			en = self.awk(path,energy)
			hess[i,i] = (2 * const) * (ep + en - (2 * exxx))
		return hess

	def awk(self,fyle,keyword):   # Works like the bash awk command
		for line in reversed(open(fyle).readlines()):
			if re.search(keyword,line) != None:
				val = re.split(keyword,line)
		return float(val[-1])
	
	def make_hessFile(self,mat):   # Makes file hessian.dat
		f = open('hessian.dat','w')
		string = "" 
		for i in range(3 * self.mol.natom):
			for j in range(3 * self.mol.natom):
				string += "{: 13.7f}".format(mat[i,j])
			string += '\n'
		f.write(string)	
		f.close()

if __name__ == "__main__":
	template = open('../../extra-files/template.dat','r').read()
	hess = Hessian('../../extra-files/molecule.xyz')
	directory = 'DISPS'
	command = 'psi4'
	disp_size = 0.005
	energy_prefix = '@DF-RHF Final Energy:'
	sleep_time = 0
	hess.generate_inputs(template,disp_size,directory)
	hess.run_jobs(command,directory,sleep_time)
	newHess = hess.build_hessian(energy_prefix,disp_size,directory)
	hess.make_hessFile(newHess)
        # Call back to Project 1
	freq = Frequency('../../extra-files/hessian.dat','../../extra-files/molecule.xyz')
	freq.frequencies(freq.mol,freq.hess)
