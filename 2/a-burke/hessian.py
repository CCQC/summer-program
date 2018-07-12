#Author: Ally Burke

import re
import os
import sys
import numpy as np
sys.path.append('../../extra-files')
sys.path.append('../../0/a-burke')
sys.path.append('../../1/a-burke')
from molecule import Molecule
from project1 import frequencies
import subprocess

#Read in molecule.xyz
with open('../../extra-files/molecule.xyz') as geometry:
	geometry = geometry.read()

#Build a Molecule object from molecule.xyz
mol = Molecule(geometry, "Angstrom")
mol.to_bohr()

#Build input file template (template) from template.dat
temp = open('template.dat','r')

template = ''
template += 'memory 256 mb' + '\n'
for t in temp:
	template += t
template = template.replace('{','',1)
template = template.replace('}}','}')
copy = template
N = int(mol.natom)

#SIMPLIFY GENERATE_INPUTS FUNCTION


directory_list = []	

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):
	#Input file for reference geometry
	refgeom = "units " + mol.units + '\n' 
	refgeom += str(mol.labels[0]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(mol.geom[0][0],mol.geom[0][1],mol.geom[0][2]) + '\n'	
	refgeom += str(mol.labels[1]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(mol.geom[1][0],mol.geom[1][1],mol.geom[1][2]) + '\n'	
	refgeom += str(mol.labels[2]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(mol.geom[2][0],mol.geom[2][1],mol.geom[2][2]) 	
	dire = directory + '/'
	
	
	if not os.path.exists(dire):
		os.mkdir(dire)
	if not os.path.exists(dire + 'd/'):
		os.mkdir(dire + 'd/')
		directory_list.append('d/')
		file = open(dire + 'd/input.dat','w')
		template = template.replace('{:s}',refgeom)
		file.write(template)
		file.close()
	

	#Positive single displacements
	for n in range(3*N):
		template = copy
		tempdisp = []
		disp = []
		for i in range(3):
			for j in range(3):
				disp.append(float(mol.geom[i][j]))
		tempdisp = disp
		newgeom = ''
		tempdisp[n] += disp_size
		tempdisp = np.array(tempdisp)
		tempdisp = tempdisp.reshape((3,N))
		newgeom = ''
		newgeom += "units " + mol.units + '\n'
		newgeom += str(mol.labels[0]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[0][0],tempdisp[0][1],tempdisp[0][2]) + '\n' 
		newgeom += str(mol.labels[1]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[1][0],tempdisp[1][1],tempdisp[1][2]) + '\n'
		newgeom += str(mol.labels[2]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[2][0],tempdisp[2][1],tempdisp[2][2])
		name = dire + 'd+' + str(n) + '/'
		directory_list.append('d+' + str(n) + '/')
		if not os.path.exists(name):
			os.mkdir(name)
			file = open(name + 'input.dat','w')
			template = template.replace('{:s}',newgeom)
			file.write(template)
			file.close()

	#Negative single displacements
	for n in range(3*N):
		template = copy
		tempdisp = []
		disp = []
		for i in range(3):
			for j in range(3):
				disp.append(float(mol.geom[i][j]))
		tempdisp = disp
		newgeom = ''
		tempdisp[n] -= disp_size
		tempdisp = np.array(tempdisp)
		tempdisp = tempdisp.reshape((3,N))
		newgeom = ''
		newgeom += "units " + mol.units + '\n'
		newgeom += str(mol.labels[0]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[0][0],tempdisp[0][1],tempdisp[0][2]) + '\n' 
		newgeom += str(mol.labels[1]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[1][0],tempdisp[1][1],tempdisp[1][2]) + '\n'
		newgeom += str(mol.labels[2]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[2][0],tempdisp[2][1],tempdisp[2][2])
		name = dire + 'd-' + str(n) + '/'
		directory_list.append('d-' + str(n) + '/')
		if not os.path.exists(name):
			os.mkdir(name)
			file = open(name + 'input.dat','w')
			template = template.replace('{:s}',newgeom)
			file.write(template)
			file.close()

	#Positive double displacements
	
	for i in range(3*N):
		for j in range(i+1,3*N):
			template = copy
			tempdisp = []
			disp = []
			for l in range(3):
				for k in range(3):
					disp.append(float(mol.geom[l][k]))
			tempdisp = disp
			tempdisp[i] += disp_size
			tempdisp[j] += disp_size
			tempdisp = np.array(tempdisp)
			tempdisp = tempdisp.reshape((3,N))
			newgeom = ''
			newgeom += "units " + mol.units + '\n'
			newgeom += str(mol.labels[0]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[0][0],tempdisp[0][1],tempdisp[0][2]) + '\n' 
			newgeom += str(mol.labels[1]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[1][0],tempdisp[1][1],tempdisp[1][2]) + '\n'
			newgeom += str(mol.labels[2]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[2][0],tempdisp[2][1],tempdisp[2][2])
			name = dire + 'd+' + str(i) + str(j) + '/'
			directory_list.append('d+' + str(i) + str(j) + '/')
			if not os.path.exists(name):
				os.mkdir(name)
				file = open(name + 'input.dat','w')
				template = template.replace('{:s}',newgeom)
				file.write(template)
				file.close()

				
	#Negative double displacements

	for i in range(3*N):
		for j in range(i+1,3*N):
			template = copy
			tempdisp = []
			disp = []
			for l in range(3):
				for k in range(3):
					disp.append(float(mol.geom[l][k]))
			tempdisp = disp
			tempdisp[i] -= disp_size
			tempdisp[j] -= disp_size
			tempdisp = np.array(tempdisp)
			tempdisp = tempdisp.reshape((3,N))
			newgeom = ''
			newgeom += "units " + mol.units + '\n'
			newgeom += str(mol.labels[0]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[0][0],tempdisp[0][1],tempdisp[0][2]) + '\n' 
			newgeom += str(mol.labels[1]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[1][0],tempdisp[1][1],tempdisp[1][2]) + '\n'
			newgeom += str(mol.labels[2]) + '{:>10.5f}{:>10.5f}{:>10.5f}'.format(tempdisp[2][0],tempdisp[2][1],tempdisp[2][2])
			name = dire + 'd-' + str(i) + str(j) + '/'
			directory_list.append('d-' + str(i) + str(j) + '/')
			if not os.path.exists(name):
				os.mkdir(name)
				file = open(name + 'input.dat','w')
				template = template.replace('{:s}',newgeom)
				file.write(template)
				file.close()

def run_jobs(mol, command = 'psi4', directory = "DISPS"):
	name = '/opt/vulcan/bin/vulcan  submit debug.q psi4@1.0'

	for d_file in directory_list:	
		dirname = '/home/vulcan/ab99072/summer-program/2/a-burke/DISPS/' + d_file
		os.chdir(dirname)
		os.system(name)
		os.chdir('..')
	
def grab_energy(mol, energy_prefix = '@DF-RHF Final Energy:'):
	f = open('output.dat').read()
	x = f.find(energy_prefix)
	x += len(energy_prefix)
	x += 3
	z = '' 
	for i in range(x,len(f)):
		if f[i] == '\n':
			break
		else:
			z += f[i]
	z = float(z)
	return z

def build_hessian(mol, energy_prefix = '@DF-RHF Final Energy:', disp_size = 0.005, directory = "DISPS"):
	disps = '/home/vulcan/ab99072/summer-program/2/a-burke/' + directory + '/'
	os.chdir(disps + 'd/')
	d = grab_energy(mol, energy_prefix)
	hessian = np.zeros((3*N,3*N))
	for i in range(3*N):
		for j in range(3*N):
			#Complete diagonal of hessian
			if i == j:
				os.chdir(disps + 'd+' + str(i))
				dp = grab_energy(mol)
				os.chdir(disps + 'd-' + str(i))
				dm = grab_energy(mol)
				h = (dp + dm - 2*d)/(disp_size**2)
				hessian[i,j] += h
			#Complete the rest of the hessian
			else:
				k = sorted([i,j])
				os.chdir(disps + 'd+' + str(k[0]) + str(k[1]))
				dpp = grab_energy(mol)
				os.chdir(disps + 'd-' + str(k[0]) + str(k[1]))
				dmm = grab_energy(mol)
				os.chdir(disps + 'd+' + str(i))
				dpn = grab_energy(mol)
				os.chdir(disps + 'd-' + str(i))
				dmn= grab_energy(mol)
				os.chdir(disps + 'd+' + str(j))
				dnp = grab_energy(mol)
				os.chdir(disps + 'd-' + str(j))
				dnm = grab_energy(mol)
				h = (dpp + dmm - dpn -dmn - dnp - dnm + 2*d)/(2*disp_size**2)
				hessian[i,j] += h

	os.chdir('..')	
	np.savetxt('hessian.dat', hessian,"%15.7f"," ","\n")
		
#import your frequencies from Project 1 and use it to calculate frequencies and normal modes from the return value of build_hessian
	
	freq = frequencies(mol,hessian)
	
	return hessian

generate_inputs(mol,template,disp_size = 0.005,directory = 'DISPS')
run_jobs(mol, command = 'psi4', directory = "DISPS")
build_hessian(mol, energy_prefix = '@DF-RHF Final Energy:', disp_size = 0.005, directory = "DISPS")
