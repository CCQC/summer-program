#!/usr/bin/python

import os
import re
import sys
sys.path.insert(0, '../../0/aewiens/')

from molecule import Molecule
import numpy as np

h = 0.005

####  Read in geometry  ####

f = open("../../1/extra-files/molecule.xyz","r").readlines()
mol = Molecule(f)
mol.bohr()
N = mol.__len__()
coords = []
for i in range(N):
	for j in range(N):
		coords.append(mol.geom[i][j])


####  function: make and run input files  #### 

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

#### function: find energy from output files ####

def E(i,j,hi,hj):
   dirname = "X%dX%d_%d%d" % (i,j,hi,hj)
   f = open("%s/output.dat" % dirname, "r").readlines()
   lines = (' ').join(f)
   line = re.findall("Total Energy\s=\s+-\d+.\d+",lines)
   return float(line[0].split()[3])

####  Run reference configuration    ####

#get_energy("X0X0_00",mol.atoms,coords)
E0 = E(0,0,0,0)

####   Run single displacements   ####

for i in range(3*N):
   forward = "X%dX0_10" % i
   coords_f = [j for j in mol.coords]
   coords_f[i] +=  h
   #get_energy(forward,mol.atoms,coords_f)

   reverse = "X%dX0_-10" % i
   coords_r = [j for j in mol.coords]
   coords_r[i] -=  h
   #get_energy(reverse,mol.atoms,coords_r)


####   Run double displacements    ######

for i in range(3*N):
   for j in range(i):
      f = "X%dX%d_11" % (i,j)
      coords_f2 = [k for k in mol.coords]
      coords_f2[i] += h
      coords_f2[j] += h 
      #get_energy(f,mol.atoms,coords_f2)

      r = "X%dX%d_-1-1" % (i,j)
      coords_r2 = [l for l in mol.coords]
      coords_r2[i] -= h
      coords_r2[j] -= h 
      #get_energy(r,mol.atoms,coords_r2)

H = np.zeros((3*N,3*N))
for i in range(3*N):
   a = E(i,0,1,0)
   b = E(i,0,-1,0)
   c = 2*E0
   print ((a + b)-c)/(h**2)
#   print (a + b - c)/(h**2)
   #H[i,i]= (E(i,0,1,0)+E(i,0,-1,0)-2*E0)/(h**2)
   #for j in range(0,i):
   #   H[i,j] = H[j,i] = (E(i,j,1,1)+E(i,j,-1,-1)-E(i,0,1,0)-E(j,0,1,0)-E(j,0,-1,0)-E(i,0,-1,0)+2*E0)/(h**2)

#np.savetxt("hessian.txt",H,"%15.7f"," ","\n")
