#part 2 - use this function to create the hessian from your data

import sys
sys.path.append("../../0/rkm79791")
sys.path.append("../../1/rkm79791")
from p0 import Molecule
import numpy as np
import os
import shutil
import subprocess
import re
from project1 import freq


##################
directory = "DISPS"
disp_size = 0.005
with open('molecule.xyz', 'r') as f: 
    f = f.readlines()
    mole = [x.replace('\n', '') for x in f]    


mol = Molecule("molecule.xyz")
mol.to_bohr()
N, geom, atoms = int(mol.num_atoms), mol.geom, mol.atom_labels
geom = np.array(geom)

##################


def energy(mol):
    with open("output.dat") as info:
        lines = info.readlines()
        for p in lines:
            if "@DF-RHF Final Energy:" in p:
                q = p.split()
                E = float(q[-1]) 
                return E

                
def make_hessian(mol, directory, disp_size):
    os.chdir(str(os.getcwd()) + "/" + directory + "/" + "DISP0")
    N = int(mol.num_atoms)
    E1 = energy(mol)
    os.chdir('..')
    Hess = np.zeros((3*N, 3*N))
    for m in range(3*N):
        for n in range(3*N):  #3N x 3N matrix....
            #we can get the "diagonal" of the matrix when i = j, so let's do that
            if m == n:      
                dire1 = str(os.getcwd()) + "/" + directory + str("+") + str(m) + "/"   
                os.chdir(dire1)
                Eposm = energy(mol)
                os.chdir('..')
                dire2 = str(os.getcwd()) + "/" + directory + str("-") + str(m) + "/"   
                os.chdir(dire2)
                Enegm = energy(mol)
                os.chdir('..')
                Eh = (Eposm + Enegm - 2*E1)/(disp_size**2)
                #Equation 2
                Hess[m,n] += Eh
            else:
                #(If m doesn't equal n, use Equation 3!!!)
                o = sorted([m,n])
                os.chdir(str(os.getcwd()) + "/" +  "DISPS" + "+" + str(o[0]) + str(o[1]))
                E1E = energy(mol)
                os.chdir("..")    
                os.chdir(str(os.getcwd()) + "/" + "DISPS" + "-" + str(o[0]) + str(o[1]))
                E2E = energy(mol)
                os.chdir("..") 
                os.chdir(str(os.getcwd()) + "/" + "DISPS" + "+" + str(m))
                E3E = energy(mol)
                os.chdir("..")
                os.chdir(str(os.getcwd()) + "/" + "DISPS" + "-" +  str(m))
                E4E = energy(mol)
                os.chdir("..")
                os.chdir(str(os.getcwd()) + "/" + "DISPS" + "+" + str(n))
                E5E = energy(mol)
                os.chdir("..")
                os.chdir(str(os.getcwd()) + "/" + "DISPS" + "-" + str(n))
                E6E = energy(mol) 
                os.chdir("..")
                                     
                H = (E1E + E2E - E3E - E4E - E5E - E6E + 2*E1)/(2*disp_size**2)
                Hess[m,n] += H
    os.chdir("..")
    np.savetxt("hessian.dat", Hess, "%18.10f", "   ", "\n")
    freq(mol)





        

gen_inputs(mol, template, disp_size = 0.005, directory = "DISPS")
run_jobs(mol, command = "psi4@1.0", directory = 'DISPS')
make_hessian(mol, directory, disp_size)
