#The second (but actually third) programming project

import sys
import os
import shutil
import errno
import math
import itertools
import numpy as np
sys.path.insert(0, '../../0/juliaanoonan')
from juliaanoonan0 import Molecule
sys.path.insert(0, '../../1/juliaanoonan')
from juliaanoonan1 import frequencies
sys.path.insert(0, '../../extra-files')

#creates a new input to add
def make_new(coords,mol):
            new = ""
            for i in range(len(mol.labels)):
                new += "{:} {:<20.10f}{:<20.10f}{:<20.10f}".format(mol.labels[i],coords[i][0],coords[i][1],coords[i][2])
                new += '\n'
            return new

def generate_inputs(mol, template, disp_size, directory):
    d = directory
    try:
        os.mkdir(d)
    except OSError as exception:
    #if it exists, so what, (for us) overwrite it
        if exception.errno != errno.EEXIST:
            raise
    #initializing things
    temp = template
    #no displacement
    if not os.path.exists("{:}/undisplaced".format(d)):
        os.mkdir("{:}/undisplaced".format(d))
    else:
        shutil.rmtree("{:}/undisplaced".format(d))
        os.mkdir("{:}/undisplaced".format(d))
    with open("{:}/undisplaced/input.dat".format(d),'w+') as f:
        f.write(temp.format(make_new(mol.geoms,mol)))
    #single displacements
    copy = np.ndarray.copy(mol.geoms)
    for i in range(3*mol.natom): #per atom x, y, z
        coords = np.ndarray.copy(copy)
        #i//3 the atom and i%3 the coord
        coords[i//3,i%3] += disp_size
        path = "{:}/{:d}{:}{:}{:d}".format(d,i,"X","+",0) #coordposition,coordposition,disp,disp
        #make directory
        if not os.path.exists(path):
            os.mkdir(path)
        else:
            shutil.rmtree(path)
            os.mkdir(path)
        #write the file
        with open("{:}/input.dat".format(path), 'w+') as g:
            g.write(temp.format(make_new(coords,mol)))
        coords[i//3,i%3] += -2*disp_size  #negative
        path = "{:}/{:d}{:}{:}{:d}".format(d,i,"X","-",0)##coordposition,coordposition,disp,disp
        if not os.path.exists(path):  #make the directory
            os.mkdir(path)
        else:
            shutil.rmtree(path)
            os.mkdir(path)
        with open("{:}/input.dat".format(path), 'w+') as g:  #write the file
            g.write(temp.format(make_new(coords,mol)))
    #double displacements
    for i in range(3*mol.natom):  #per atom x, y, z
        for j in range(i):   #always one less! matches each with every slot below
            coords = np.ndarray.copy(copy)
            coords[i//3,i%3] += disp_size   #two pos displacements
            coords[j//3,j%3] += disp_size
            path = "{:}/{:d}{:d}{:}{:}".format(d,i,j,"+","+")
            if not os.path.exists(path):  #just overwrites if exists
                os.mkdir(path)
            else:
                shutil.rmtree(path)
                os.mkdir(path)
            with open("{:}/input.dat".format(path), 'w+') as f:
                f.write(temp.format(make_new(coords,mol)))
            coords[i//3,i%3] += -2*disp_size
            coords[j//3,j%3] += -2*disp_size
            path = "{:}/{:d}{:d}{:}{:}".format(d,i,j,"-","-")
            if not os.path.exists(path):
                os.mkdir(path)
            else:
                shutil.rmtree(path)
                os.mkdir(path)
            with open("{:}/input.dat".format(path), 'w+') as v:
                v.write(temp.format(make_new(coords,mol)))

#4 walk through directories created and execute
def run_jobs(mol, command, directory):
#    os.chdir(directory)
    subdirs = [x[1] for x in os.walk(directory)]
    subdirs = subdirs[0]
    for name in subdirs:
        os.chdir("{:}/{:}".format(directory, name))
        os.system(command)
        os.chdir('../..')

#5
def yankE(path, energy_prefix="@DF-RHF Final Energy:"):
    with open("{:}/output.dat".format(path), 'r') as f:  #open the file
        search = f.readlines()  #copy it to search
    look = energy_prefix.split()  #finding final word in energy_prefix to index when we split the line
    ind = look[-1]
    for line in search:
        if energy_prefix in line:
            line = line.split()
            return float(line[line.index(ind)+1])
    print("Error! Could not find energy.")

def build_hessian(mol, energy_prefix, disp_size, directory):
    H = []
    d = directory
    for i in range(9): #9 rows
        myrow = []
        for j in range(9): #9 columns
            if i == j:  #eq 2
               el = yankE("{:}/{:d}X+0".format(d,i))+yankE("{:}/{:d}X-0".format(d,i))-2*yankE("{:}/undisplaced".format(d))
               el = el/(disp_size**2)
               myrow += [el]
            else:  #will be eq 3
                if i > j:  #the way I set up the for loops, the first num needs to be bigger, hence the ifs
                    el = yankE("{:}/{:d}{:d}++".format(d,i,j))+yankE("{:}/{:d}{:d}--".format(d,i,j))-yankE("{:}/{:d}X+0".format(d,i))-yankE("{:}/{:d}X-0".format(d,i))-yankE("{:}/{:d}X+0".format(d,j))-yankE("{:}/{:d}X-0".format(d,j))+2*yankE("{:}/undisplaced".format(d))
                    el = el/(2*(disp_size**2))
                    myrow += [el]
                if j > i:
                    el = yankE("{:}/{:d}{:d}++".format(d,j,i))+yankE("{:}/{:d}{:d}--".format(d,j,i))-yankE("{:}/{:d}X+0".format(d,i))-yankE("{:}/{:d}X-0".format(d,i))-yankE("{:}/{:d}X+0".format(d,j))-yankE("{:}/{:d}X-0".format(d,j))+2*yankE("{:}/undisplaced".format(d))
                    el = el/(2*(disp_size**2))
                    myrow += [el]
        H.append(myrow)
    H = np.array(H)
    np.savetxt("myhessian.dat",H)
    return H

#1 build molecule object
mol = Molecule('../../extra-files/molecule.xyz')

#2 build input fild template
template = open('../../extra-files/template.dat','r').read()

#6 calculate frequencies and normal modes from my Hessian
if __name__ == "__main__":
    generate_inputs(mol, template, disp_size = (0.005*0.529177249), directory = "DISPS")
    run_jobs(mol, command = "psi4", directory = "DISPS")
    x = build_hessian(mol, "@DF-RHF Final Energy:", disp_size = (0.005), directory = "DISPS")
    frequencies(mol,x)

