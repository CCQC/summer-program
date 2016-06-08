#! /usr/bin/python3

import os
import numpy as np
from molecule import Molecule

def mkinput(path,inp):
    if not os.path.exists(path):
        os.makedirs(path)
    new = open(path + '/' + 'input.dat', 'w')
    new.write(inp)
    new.close()

def grab_energy(energy_prefix, output):
    with open(output) as out:
        for line in out:
            if energy_prefix in line:
                r = line.replace(energy_prefix, '')
                break
    return float(r)

def write(matrix):
    println = '' 
    for line in matrix: 
        for item in line:
            println += '{:>13.8f}'.format(item)
        println += '\n'
    return println

def generate_inputs(molecule, template, disp_size = 0.005, directory = 'DISPS' ):

# Reads template and store it on temp

    with open(template) as t:
        temp = ''
        for lines in t:
            temp += lines

# Creates the directory and input file for the reference structure

    path = directory
    inp = temp.format(molecule.stru())
    mkinput(path + '/ref', inp)

# Creates the directories for single displacements
    
    path_sp = path +  '/c{:d}p'
    path_sn = path +  '/c{:d}n'
    for i in range(len(molecule)):
        for j in range(3):
            n = 3*i + j
            save = molecule.geom[i][j]                        # Creates positive displacements 
            molecule.geom[i][j] += disp_size
            inp = temp.format(molecule.stru())
            mkinput(path_sp.format(n), inp)
            molecule.geom[i][j] -= 2*disp_size                # Creates negative displacements 
            inp = temp.format(molecule.stru())
            mkinput(path_sn.format(n), inp)
            molecule.geom[i][j] = save
            
# Creates the directories for double displacements

    path_dp = path + '/c{:d}c{:d}p'
    path_dn = path + '/c{:d}c{:d}n'
    for a in range(3*len(molecule)):
        b = 3*len(molecule) - 1
        a_j = int(a%3)
        a_i = int((a - a_j)/3)
        save_a = molecule.geom[a_i][a_j]
        while a < b:
            b_j = int(b%3)
            b_i = int((b - b_j)/3)
            save_b = molecule.geom[b_i][b_j]
            molecule.geom[a_i][a_j] += disp_size              # Creates positive double displacements
            molecule.geom[b_i][b_j] += disp_size
            inp = temp.format(molecule.stru()) 
            mkinput(path_dp.format(a, b), inp) 
            molecule.geom[a_i][a_j] -= 2* disp_size           # Creates negative displacements
            molecule.geom[b_i][b_j] -= 2* disp_size
            inp = temp.format(molecule.stru()) 
            mkinput(path_dn.format(a, b), inp) 
            molecule.geom[a_i][a_j] = save_a
            molecule.geom[b_i][b_j] = save_b
            b -= 1
        
def run_jobs(molecule, command = 'psi4', directory = 'DISPS'):
        
    os.chdir(directory)

# Runs the ref structure
    
    os.chdir('ref')
    os.system(command)
    os.chdir('../')
    
# Runs the single displacements

    for i in range(3*len(molecule)):
        os.chdir('c{:d}p'.format(i))
        os.system(command)
        os.chdir('../c{:d}n'.format(i))
        os.system(command)
        os.chdir('../')

# Runs the double displacements

    for a in range(3*len(molecule)):
        b = 3*len(molecule) - 1
        while a < b:
            os.chdir('c{:d}c{:d}p'.format(a, b))
            os.system(command)
            os.chdir('../c{:d}c{:d}n'.format(a, b))
            os.system(command)
            os.chdir('../')
            b -= 1

    os.chdir('../')

def build_hessian(molecule, e = '@DF-RHF Final Energy:', disp_size = 0.005, directory = 'DISPS'):
    h = disp_size
    hessian = np.zeros((3*len(molecule), 3*len(molecule)))
    E0 = grab_energy(e, directory + '/ref/output.dat')
    for x in range(len(hessian)):
        for y in range(len(hessian))[x:]:
            if x == y:                                                                        # Second derivative with recpect to the same variable
                Ep = grab_energy(e, directory + '/c{:d}p/output.dat'.format(x))
                En = grab_energy(e, directory + '/c{:d}n/output.dat'.format(x))
                hessian[x][y] = (Ep + En - 2*E0)/(h**2) 
            else:                                                                             # Second derivative with respecto to two different variables
                Epx = grab_energy(e, directory + '/c{:d}p/output.dat'.format(x))
                Epy = grab_energy(e, directory + '/c{:d}p/output.dat'.format(y))
                Enx = grab_energy(e, directory + '/c{:d}n/output.dat'.format(x))
                Eny = grab_energy(e, directory + '/c{:d}n/output.dat'.format(y))
                Epp = grab_energy(e, directory + '/c{:d}c{:d}p/output.dat'.format(x, y))
                Enn = grab_energy(e, directory + '/c{:d}c{:d}n/output.dat'.format(x, y))
                hessian[x][y] = (Epp + Enn - Epx - Enx - Epy - Eny + 2*E0) / (2*h**2)
                hessian[y][x] = hessian[x][y]
    out = open('hessian.dat', 'w')
    out.write(write(hessian))
    out.close()
    
water = Molecule('water.xyz')
water.bohr()
generate_inputs(water, 'template.dat')
run_jobs(water)
build_hessian(water)
    
