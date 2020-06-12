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




directory = "DISPS"
disp_size = 0.005
with open('molecule.xyz', 'r') as f: 
    f = f.readlines()
    mole = [x.replace('\n', '') for x in f]    


mol = Molecule("molecule.xyz")
mol.to_bohr()
N, geom, atoms = int(mol.num_atoms), mol.geom, mol.atom_labels
geom = np.array(geom)

#Make geom template like in project 1:

def make_temp(geom, atoms):
        addon = "units bohr\n"
        for i in range(len(atoms)):
            x,y,z = geom[i]
            addon += "{:} {:>15.10f}{:>15.10f}{:>15.10f}".format(atoms[i], x, y, z)
            if i < (len(atoms)-1):
                addon += "\n"
        return addon

######################

template = ""
with open("template.dat", "r") as file:
    temp = file.readlines()
    for i in temp:
        template += i

#######################

d_list = []


def gen_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):
    
    dire = directory 
    
    if not os.path.exists(dire):
        os.mkdir(dire)
    if not os.path.exists(dire + "/DISP0"):
        os.mkdir(dire+"/DISP0")
    d_list.append("/DISP0")
    #Make initial dirs for zeroth displacements
    
    os.chdir(dire +"/DISP0")
    ######
    
    if not os.path.exists("input.dat"):
        temp2 = template.replace("{s}", make_temp(geom,atoms))
        #print(template)
        with open("input.dat", "w") as file: 
            #template = template.replace("{:s}",make_temp(geom,atoms))
            file.writelines(temp2)
            file.close 
    os.chdir("..")
    
    ######
    copy = np.ndarray.copy(mol.geom)


#single displacements    
    for i in range(3*N):
            coords = np.ndarray.copy(copy)
            coords[i//3,i%3] += disp_size
            #make directory
            dire1 = str("+") + str(i) + "/"    #keep track - this is positive displacement
            dire2 = directory + dire1
            d_list.append(dire2)
            disp_ary1 = np.array(coords).reshape(3,N)
            #print(disp_ary)
            if not os.path.exists(dire2):
                os.mkdir(dire2)
                os.chdir(dire2)
                newgeo = "units bohr \n"
                for j in range(N):
                    a,b,c = disp_ary1[j]
                    #print(a)
                    newgeo += "{:} {:>15.10f}{:>15.10f}{:>15.10f}".format(atoms[j], a, b, c)
                    if j < (len(atoms)-1):
                        newgeo += "\n"
                #print(newgeo)
                temp1 = template.replace("{s}", newgeo)                  
                #print(temp1)
                with open("input.dat", "w") as file: 
                    file.writelines(temp1)
                    file.close
                os.chdir("..")
    
    
    for i in range(3*N):
            coords = np.ndarray.copy(copy)
            coords[i//3,i%3] -= disp_size
            #make directory
            dire1 = str("-") + str(i) + "/"    #keep track - this is negative displacement
            dire2 = directory + dire1
            d_list.append(dire2)
            disp_ary2 = np.array(coords).reshape(3,N)
            if not os.path.exists(dire2):
                os.mkdir(dire2)
                os.chdir(dire2)
                newgeo = "units bohr\n"
                for j in range(N):
                    a,b,c = disp_ary2[j]
                    newgeo += "{:} {:>15.10f}{:>15.10f}{:>15.10f}".format(atoms[j], a, b, c)
                    if j < (len(atoms)-1):
                        newgeo += "\n"
                temp1 = template.replace("{s}", newgeo)                  
                with open("input.dat", "w") as file: 
                    file.writelines(temp1)
                    file.close
                os.chdir("..")
    
    
    #######Try some double displacements
    
    
    for i in range(3*N):
        for j in range(i+1,3*N):
                coords = np.ndarray.copy(copy)
                coords[i//3,i%3] += disp_size
                coords[j//3,j%3] += disp_size
                #make directory
                dire1 = str("+") + str(i) + str(j) + "/"    #keep track - this is positive displacement
                dire2 = directory + dire1
                d_list.append(dire2)
                disp_ary3 = np.array(coords).reshape(3,N)
                if not os.path.exists(dire2):
                    os.mkdir(dire2)
                    os.chdir(dire2)
                    newgeo = "units bohr\n"
                    for r in range(N):
                        a,b,c = disp_ary3[r]
                        newgeo += "{:} {:>15.10f}{:>15.10f}{:>15.10f}".format(atoms[r], a, b, c)
                        if r < (len(atoms)-1):
                            newgeo += "\n"
                    temp1 = template.replace("{s}", newgeo)                  
                    with open("input.dat", "w") as file: 
                        file.writelines(temp1)
                        file.close
                    os.chdir("..")
            
    
    for i in range(3*N):
        for j in range(i+1,3*N):
                coords = np.ndarray.copy(copy)
                coords[i//3,i%3] -= disp_size
                coords[j//3,j%3] -= disp_size
                #make directory
                dire1 = str("-") + str(i) + str(j) + "/"    #keep track - this is negative displacement
                dire2 = directory + dire1
                d_list.append(dire2)
                disp_ary4 = np.array(coords).reshape(3,N)
                if not os.path.exists(dire2):
                    os.mkdir(dire2)
                    os.chdir(dire2)
                    newgeo = "units bohr\n"
                    for r in range(N):
                        a,b,c = disp_ary4[r]
                        newgeo += "{:} {:>15.10f}{:>15.10f}{:>15.10f}".format(atoms[r], a, b, c)
                        if r < (len(atoms)-1):
                            newgeo += "\n"
                    temp1 = template.replace("{s}", newgeo)                  
                    with open("input.dat", "w") as file: 
                        file.writelines(temp1)
                        file.close
                    os.chdir("..")


#run jobs  


def run_jobs(mol, command = "psi4@1.0", directory = 'DISPS'):
    loc = str(os.getcwd()) + "/"
    c = '/opt/vulcan/bin/vulcan  submit gen3.q ' + command
    os.chdir(loc)
    for i in d_list:
        d = loc + i   #directory
        os.chdir(d)   
        os.system(c)  #execute command to gen3
        os.chdir('..')  

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
