#! /usr/bin/python3

import os
import numpy as np
import sys
sys.path.insert(0, '../../0/aroeira')
sys.path.insert(0, '../../1/aroeira')
sys.path.insert(0, '../../1/extra-files')
from molecule import Molecule
from frequencies import Hessian

def mkinput(path,inp):                                           # Creates the path and inside it creates a input file
    if not os.path.exists(path):
        os.makedirs(path)
    new = open(path + '/' + 'input.dat', 'w')
    new.write(inp)
    new.close()

def grab_energy(energy_prefix, output):                          # Extracts the energy value from a psi4 output
    with open(output) as out:
        for line in out:
            if energy_prefix in line:
                r = line.replace(energy_prefix, '')
                break
    return float(r)

def mwrite(matrix):                                              # Write a matrix as a formated string
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
    for i in range(len(molecule)):                            # Iterate through the number of atoms N
        for j in range(3):                                    # For each atom iterate through the three coordinates
            n = 3*i + j                                       # Determine the absolute index of a given coordinate in order to organize the directories
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
    for a in range(3*len(molecule)):                          # Iterate through the number of coordinates using absolute index
        b = 3*len(molecule) - 1                               # Start b as the last coordinate
        a_j = int(a%3)                                        # Determine the geometry index of the absolute index 'a'
        a_i = int((a - a_j)/3)
        save_a = molecule.geom[a_i][a_j]                      # Store the coordinate 'a' of the reference structure
        while a < b:                                          # Iterate through all unique pairs
            b_j = int(b%3)                                    # Determine the geometry index of the absolute index 'b'
            b_i = int((b - b_j)/3)
            save_b = molecule.geom[b_i][b_j]                  # Store the coordinate 'b' of the reference strucutre
            molecule.geom[a_i][a_j] += disp_size              # Creates positive double displacements
            molecule.geom[b_i][b_j] += disp_size
            inp = temp.format(molecule.stru())                # Makes the input 
            mkinput(path_dp.format(a, b), inp) 
            molecule.geom[a_i][a_j] -= 2* disp_size           # Creates negative displacements
            molecule.geom[b_i][b_j] -= 2* disp_size
            inp = temp.format(molecule.stru()) 
            mkinput(path_dn.format(a, b), inp)                # Makes the input
            molecule.geom[a_i][a_j] = save_a                  # Restore the alterated geometry to their original value
            molecule.geom[b_i][b_j] = save_b
            b -= 1                                            # Update b
        
def run_jobs(molecule, command = 'psi4', directory = 'DISPS'):
        
    os.chdir(directory)

# Runs the ref structure
    
    os.chdir('ref')
    os.system(command)
    os.chdir('../')
    
# Runs the single displacements

    for i in range(3*len(molecule)):                        # Iterates through the number of coordinates
        os.chdir('c{:d}p'.format(i))
        os.system(command)
        os.chdir('../c{:d}n'.format(i))
        os.system(command)
        os.chdir('../')

# Runs the double displacements

    for a in range(3*len(molecule)):                        
        b = 3*len(molecule) - 1
        while a < b:                                        # Iterates through all unique pair of coordinates, except a = b
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
            if x == y:                                                                        # Second derivative with respect to the same variable
                Ep = grab_energy(e, directory + '/c{:d}p/output.dat'.format(x))
                En = grab_energy(e, directory + '/c{:d}n/output.dat'.format(x))
                hessian[x][y] = (Ep + En - 2*E0)/(h**2) 
            else:                                                                             # Second derivative with respect to two different variables
                Epx = grab_energy(e, directory + '/c{:d}p/output.dat'.format(x))
                Epy = grab_energy(e, directory + '/c{:d}p/output.dat'.format(y))
                Enx = grab_energy(e, directory + '/c{:d}n/output.dat'.format(x))
                Eny = grab_energy(e, directory + '/c{:d}n/output.dat'.format(y))
                Epp = grab_energy(e, directory + '/c{:d}c{:d}p/output.dat'.format(x, y))
                Enn = grab_energy(e, directory + '/c{:d}c{:d}n/output.dat'.format(x, y))
                hessian[x][y] = (Epp + Enn - Epx - Enx - Epy - Eny + 2*E0) / (2*h**2)
                hessian[y][x] = hessian[x][y]
    out = open('hessian.dat', 'w')
    out.write(mwrite(hessian))
    out.close()
    
water = Molecule('../../extra-files/molecule.xyz')
water.bohr()
generate_inputs(water, '../../extra-files/template.dat')
run_jobs(water)
build_hessian(water)
H = Hessian(water, 'hessian.dat')
H.output()

