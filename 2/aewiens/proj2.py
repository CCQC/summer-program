#!/usr/bin/python

import os
import sys

sys.path.insert(0, '../../0/aewiens/')
from molecule import Molecule

############################
####  Read in geometry  ####
############################

f = open("../../1/extra-files/molecule.xyz","r").readlines()
mol = Molecule(f)
#mol.bohr()
geom = mol.geom
atoms = mol.atoms
N = mol.__len__()
coords = []
for i in range(N):
	for j in range(N):
		coords.append(geom[i][j])


############################
####  Make input files  #### 
############################

def get_energy(dirname,labels,coords):
	os.mkdir(dirname)
	f = open("%s/input.dat" % dirname,"w")
	f.write("set basis cc-pVDZ\n\nmolecule h20 {\n0 1\n")
	for i in range(N):
		f.write("%s %f %f %f\n"  % (labels[i], coords[3*i], coords[3*i+1], coords[3*i+2]))
	f.write("units bohr\n}\n\nenergy('scf')")
	f.close()
	os.chdir(dirname)
	os.system('psi4')
	os.chdir('..')

######################################
####  Reference Configuration    ####
######################################

get_energy("X0X0_00",atoms,coords)


##################################
####   Single Displacements   ####
##################################

for i in range(3*N):
	forward = "X%dX0_10" % i
	coords_f = coords
	coords_f[i] = coords_f[i] + 0.005
	get_energy(forward,atoms,coords_f)


	reverse = "X%dX0_-10" % i
	coords_r = coords
	coords_r[i] = coords_r[i] - 0.005
	get_energy(reverse,atoms,coords_r)


#####################################
####   Double Displacements    ######
#####################################

###########   BUG   #################

for i in range(3*N):
	for j in range(i+1,3*N):
		forward2 = "X%dX%d_11" % (i,j) 
		coords_f2 = coords
		coords_f2[i] = coords_f2[i] + 0.005
		coords_f2[j] = coords_f2[i] + 0.005
		get_energy(forward2,atoms,coords_f2)
"""
		reverse2 = "X%dX%d_-1-1" % (i,j) 
		coords_r2 = coords
		coords_r2[i] = coords_r2[i] - 0.005
		coords_r2[j] = coords_r2[j] - 0.005
		#get_energy(reverse2,atoms,coords_r2)
"""
		
