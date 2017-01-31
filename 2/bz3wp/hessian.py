# Computes Hessian matrix from single-point energies at displaced geometries by finite differences.

import sys 
sys.path.insert(0, '../../0/bz3wp')
sys.path.insert(0, '../../1/bz3wp')
from molecule import Molecule
import frequencies
import numpy as np
import os, subprocess
from glob import glob


def generate_inputs(mol, template, disp_size = 0.005, directory = 'DISPS'):         # generate inputs with all possible displacements

    mol.to_bohr()
    natom, geom, labels = mol.natom, mol.geom, mol.labels
        
    def write_input(coords, coordd, ks, kd, disp_geom):                             # function to format how to write input
        path = directory + '/{:d}{:d}_{:d}{:d}'.format(coords, coordd, ks, kd) 
        os.makedirs(path, exist_ok=True) 
        with open(path + '/input.dat', 'w') as f:
            f.write(template.format(str(disp_geom) + '\n units bohr'))
   
    def geom_conv(labels, geom):                                                   # function to convert geometry into xyz format
        geom_mrx = ''
        for atom, line in zip(labels, geom):
            x,y,z=line
            geom_mrx += "{:2s}{:>10.5f}{:>10.5f}{:>10.5f}\n".format(atom,x,y,z)
        return geom_mrx
    
    geom_form = geom_conv(labels, geom)
    write_input(0, 0, 0, 0, geom_form)                                              # writing input for initial geometry
    
    # writing input for all single displacements
    for coords in range(3*natom):
        for k in [-1, 1]:
            disps = np.zeros(3*natom)
            disps[coords] += k*disp_size
            disps = disps.reshape(geom.shape)
            disps_geom = disps + geom
            disp_geom = geom_conv(labels, disps_geom)   
            write_input(coords, 0, k, 0, disp_geom)
      
            # writing input for double displacements
            for coordd in range(coords):                                            # only for displacements in lower triangle of matrix since will be symmetric
                dispd = np.zeros(3*natom)
                dispd[coordd] += k*disp_size
                dispd = dispd.reshape(disps_geom.shape)
                dispd_geom = disps_geom + dispd  
                disp_geomd = geom_conv(labels, dispd_geom)
                write_input(coords, coordd, k, k, disp_geomd)
        
def run_jobs(mol, command = 'psi4', directory = 'DISPS'):                           # goes through all inputs and runs with psi4

    for direc in glob('DISPS/*'):
        subprocess.check_call(command + ' -i ' + direc + '/input.dat', shell=True)

def build_hessian(mol, ep, disp_size = 0.005,  directory = 'DISPS'):                # builds hessian from single point energies
    natom = mol.natom
    def E(coords, coordd, ks, kd, ep):                                              # gets energies for different displacements
        if coordd > coords and ks != 0 and kd !=0:                                  # assigning energies to correct displacements
            coords, coordd, ks, kd = coordd, coords, ks, kd
        if coordd > coords and kd == 0: 
            coords, coordd, ks, kd = coords, 0, ks, 0
        if coordd > coords and ks ==0:
            coords, coordd, ks, kd = coordd, 0, kd, 0 
        if coords > coordd and ks == 0:
            coords, coordd, ks, kd = coordd, 0, kd, 0
        if coords > coordd and kd ==0:
            coords, coordd, ks, kd = coords, 0, ks, 0
        path = directory + '/{:d}{:d}_{:d}{:d}'.format(coords, coordd, ks, kd)
        lines = open(path + '/input.out').readlines()
        energy = 0.0
        for line in reversed(lines):
            if line[:23] == ep:
                energy = float(line.split()[-1])
                return energy
        raise Exception('Cannot find energy in ' + path)
    
    H = np.zeros((3*natom, 3*natom))

    for A in range(3*natom):                                                       # calculates matrix elements of hessian
        for B in range(3*natom):
            if A == B:
                H[A,A] = (E(A,0,+1,0,ep) + E(A,0,-1,0,ep) - 2 * E(0,0,0,0,ep))/(disp_size**2)
            else:
                H[A,B] = (E(A,B,+1,+1,ep) + E(A,B,-1,-1,ep) - E(A,B,+1,0,ep) - E(A,B,-1,0,ep) \
                          - E(A,B,0,+1,ep) - E(A,B,0,-1,ep) + 2*E(0,0,0,0,ep))/(2*(disp_size**2))
    np.savetxt('hessian.dat', H)                                                   # writes hessian to hessian.dat
    return H

# testing
if __name__=='__main__':

    mol = Molecule.from_file('../../extra-files/molecule.xyz')     # building molecule from molecule.xyz
    template = open('../../extra-files/template.dat').read()       # reads initial template
    generate_inputs(mol, template, 0.005, 'DISPS')
    run_jobs(mol, '/Users/boyi/bin/psi4/obj/stage/usr/local/bin/psi4', 'DISPS')
    build_hessian(mol,'  @DF-RHF Final Energy:', 0.005, 'DISPS')
    frequencies.get_freqs(mol, 'hessian.dat')                       # calling frequencies function
