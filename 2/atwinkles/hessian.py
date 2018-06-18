#!/usr/bin/python

import os
import re
import numpy as np

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):
    
    mol = mol.copy(mol)
    mol.to_bohr()

    def write_input(coord1, coord2, k1, k2, disp_mol):
        path = directory + "/x{:d}x{:d}_{:d}{:d}".format(coord1, coord2, k1, k2)
        if not os.path.isdir(path): os.makedirs(path)
        inpf = open(path + "/input.dat",'w+')
        inpf.write(template.format(str(disp_mol)))
        inpf.close()

    # initial
    write_input(0,0,0,0, mol)

    # displacements
    for coord1 in range(3*mol.natom):
        for k in [-1,+1]:
            # single displacement
            disp1  = np.zeros(3*mol.natom)
            disp1[coord1] += k * disp_size
            disp1_mol = mol.copy(mol)
            disp1_mol.displace(disp1.reshape(mol.natom,3))
            write_input(coord1, 0, k, 0, disp1_mol)

            # double displacement
            for coord2 in range(coord1):
                disp2 = np.zeros(3*mol.natom)
                disp2[coord2] += k * disp_size
                disp2_mol = disp1_mol.copy(disp1_mol)
                disp2_mol.displace(disp2.reshape(mol.natom,3))
                write_input(coord1, coord2, k, k, disp2_mol)

def run_jobs(mol, command = "psi4", directory = "DISPS"):

    def run(coord1, coord2, k1, k2):
        os.chdir(directory + "/x{:d}x{:d}_{:d}{:d}".format(coord1, coord2, k1, k2))
        print("running ... " + ("x{:d}x{:d}_{:d}{:d}").format(coord1, coord2, k1, k2))
        os.system(command)
        os.chdir('../..')

    run(0, 0, 0, 0)
    for coord1 in range(3 * mol.natom):
        for k in [-1, +1]:
            run(coord1, 0, k, 0)
            for coord2 in range(coord1):
                run(coord1, coord2, k, k)

def build_hessian(mol, energy_prefix, disp_size = 0.005, directory = "DISPS"):

    hdim = 3 * mol.natom

    def E(coord1, coord2, k1, k2):
        coord1 = coord1 if not k1 is 0 else 0
        coord2 = coord2 if not k2 is 0 else 0
        if coord2 > coord1 or k1 is 0: coord1, coord2, k1, k2 = coord2, coord1, k2, k1
        path = directory + "/x{:d}x{:d}_{:d}{:d}".format(coord1, coord2, k1, k2)
        try:
            outf = open(path + "/output.dat").read()
            for match in re.finditer(energy_prefix + '\s+(-?\d+\.\d+)', outf): pass
            energy = float(match.group(1))
            return energy
        except:
            raise Exception("Failed to parse energy in {:s}".format(path))

    H = np.zeros((hdim, hdim))

    for A in range(hdim):
        for B in range(hdim):
            if A == B:
                H[A,A] = (E(A, A, +1, 0) + E(A, A, -1, 0) - 2 * E(A, A, 0, 0))/(disp_size**2)
            else:
                H[A,B] = (E(A,B,+1,+1) + E(A,B,-1,-1) - E(A,B,+1,0) - E(A,B,-1,0) - E(A,B,0,+1) - E(A,B,0,-1) + 2 * E(A,B,0,0))/(2 * disp_size**2)

    np.savetxt("hessian.dat", H)

    return H



if __name__ == "__main__":
    import sys
    sys.path.insert(0, "../../0/atwinkles/")
    from molecule import Molecule
    
    # builds the molecule
    mol = Molecule(open("../../extra-files/molecule.xyz").read())

    # reads in the initial template
    template = open("../../extra-files/template.dat").read()

    # generate input files
    generate_inputs(mol, template)

    # runs jobs
    run_jobs(mol)

    # build hessian
    H = build_hessian(mol, '@DF-RHF Final Energy:')
