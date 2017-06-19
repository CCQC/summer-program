import sys
sys.path.append('../../0/zachglick')
sys.path.append('../../1/zachglick')

import os
import glob
import subprocess
import itertools
import numpy as np

from molecule import Molecule
import hessian_to_freqs

def generate_inputs(mol, template, disp_size = 0.005, directory = 'DISPS') :

    fpath = "%s/initial" % ( directory )
    os.makedirs(fpath, exist_ok=True)
    with open('%s/input.dat' % fpath, 'w') as f:
        f.write(template.format(get_atomstring(mol)))

    N = mol.natom
    atom_options = tuple(range(0, mol.natom))
    xyz_options = (0,1,2)
    disp_options = (-1*disp_size, disp_size)
    single_disps = itertools.product(atom_options, xyz_options, disp_options)
    double_disps = itertools.product(atom_options, xyz_options, atom_options, xyz_options, disp_options)

    for i, disp in enumerate(single_disps):
        fpath = "%s/single_%d" % ( directory, i )
        os.makedirs(fpath, exist_ok=True)
        newmol = mol.copy()
        newmol.geom[disp[0]][disp[1]] += disp[2]
        with open('%s/input.dat' % fpath, 'w') as f:
            f.write(template.format(get_atomstring(newmol)))

    for i, disp in enumerate(double_disps):
        if disp[0] > disp[2] or (disp[0] == disp[2] and disp[1] >= disp[3]) :
            continue
        fpath = "%s/double_%d" % ( directory, i )
        os.makedirs(fpath, exist_ok=True)
        newmol = mol.copy()
        newmol.geom[disp[0]][disp[1]] += disp[4]
        newmol.geom[disp[2]][disp[3]] += disp[4]
        with open('%s/input.dat' % fpath, 'w') as f:
            f.write(template.format(get_atomstring(newmol)))

def run_jobs(mol, command = 'psi4', directory = 'DISPS') :

    for fpath in glob.glob('*/*/'):
        subprocess.call(["psi4", "%s/input.dat" % fpath])

def build_hessian(mol, energy_prefix, disp_size = 0.005, directory = 'DISPS') :
    
    N = mol.natom
    atom_options = tuple(range(0, mol.natom))
    xyz_options = (0,1,2)
    disp_options = (-1*disp_size, disp_size)
    single_disps = itertools.product(atom_options, xyz_options, disp_options)
    double_disps = itertools.product(atom_options, xyz_options, atom_options, xyz_options, disp_options)
    hess = np.zeros((3*N, 3*N))
    fpath = "%s/initial/input.out" % ( directory )
    init_energy = get_energy(fpath, energy_prefix)

    for i, disp in enumerate(single_disps):
        fpath = "%s/single_%d/input.out" % ( directory, i )
        energy = get_energy(fpath, energy_prefix)
        hess_ind = 3*disp[0]+disp[1]
        
        for num in range(3*N):
            hess[hess_ind, num] -= energy
            hess[num, hess_ind] -= energy

        hess[hess_ind, hess_ind] += 3*energy

    for i, disp in enumerate(double_disps):
        if disp[0] > disp[2] or (disp[0] == disp[2] and disp[1] >= disp[3]) :
            continue
        fpath = "%s/double_%d/input.out" % ( directory, i )
        energy = get_energy(fpath, energy_prefix)
        hess_ind1 = 3*disp[0]+disp[1]
        hess_ind2 = 3*disp[2]+disp[3]

        hess[hess_ind1, hess_ind2] += energy
        hess[hess_ind2, hess_ind1] += energy

    for row in range(3*N):
        for col in range(3*N):
            if row == col:
                hess[row, col] -= 2*init_energy
                hess[row,col] /= (disp_size ** 2)
            else:
                hess[row, col] += 2*init_energy
                hess[row,col] /= (2 * (disp_size ** 2))

    with open('numerical_hessian.dat', 'w') as f:
        for i in range(3*N):
            f.write("\t".join([str(val) for val in hess[i]]) + "\n")
    return hess

def get_atomstring(mol):
    mol.to_angstrom()
    formatted_atoms = []
    for atom in mol :
        formatted_atoms.append("%s %f %f %f" % (atom[0], atom[1][0], atom[1][1], atom[1][2]))
    mol.to_bohr()
    return '\n'.join(formatted_atoms)

def get_energy(fpath, energy_prefix):
    with open(fpath) as f:
        for line in f:
            if energy_prefix in line:
                return float(line.split()[-1])
    #throw error - energy not found ?
    return 0.0

if __name__ == '__main__':

    with open('../../extra-files/molecule.xyz', 'r') as f:
        molecule = Molecule(f.read())
        molecule.to_bohr()

    with open('../../extra-files/template.dat', 'r') as f:
        template_str = f.read()

    generate_inputs(molecule, template_str)
    run_jobs(molecule)
    hessian = build_hessian(molecule, '@DF-RHF Final Energy:')
    hessian_to_freqs.solve('../../extra-files/molecule.xyz','numerical_hessian.dat')
