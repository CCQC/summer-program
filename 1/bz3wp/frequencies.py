# Reads in a Hessian matrix, computes frequencies and normal modes, and prints to a .xyz file

import sys
sys.path.insert(0, '../../0/bz3wp/')
sys.path.insert(0, '../../extra-files')
import masses
from molecule import Molecule                               # imports Molecule class from project 0
import numpy as np


def get_freqs(mol, hes = '../extra-files/hessian.dat'):

    mol.to_angstrom()                                       # converting geometry to be in angstroms
    mass = []                                               # forming a list with 1/sqrt(individual masses)
    for atom in mol.labels:                                 # for each of the x, y, x coordinates
        mass += [masses.get_mass(atom)**(-0.5)]*3

    M = np.diag(mass)                                       # putting list in diagonal matrix
   
    hessian = np.loadtxt(hes)                               # reading in hessian file

    H = M @ hessian @ M                                     # building mass weighted hessian
    
    k, q = np.linalg.eigh(H)                                #diagonalizing mass weighted hessian

    Q = M @ q                                               #un-mass-weight eigenvectors

    #determining spatial frequencies in cm^-1
    hartreetoJ = 4.3597443e-18
    amutokg = 1.6605389e-27
    bohrtom = 5.2917721e-11
    c = 2.99792458e10
    ka = k * hartreetoJ * (1/amutokg) * (1/bohrtom)**2      # ka in a.u.
    va_conv = (1/(2 * np.pi)) * (1/c)
    va = ka * va_conv                                       # va is frequencies in cm^-1
      
    # formatting for output: geometries, frequencies, and normal modes
    out = ''
    line_form = '{:2s}' + '{: >15.10f}'*6 + '\n'
    for a , k_a  in enumerate(ka): 
        if k_a < 0: 
            out += '{}\n{: >7.2f}i cm^-1\n'.format(mol.natom, np.sqrt(np.absolute(k_a))*va_conv)
        else:
            out += '{}\n{: >7.2f} cm^-1\n'.format(mol.natom, np.sqrt(k_a)*va_conv)
        for i in range(mol.natom):
            atom = mol.labels[i]
            x, y, z = mol.geom[i]
            dx, dy, dz = Q[3*i: 3*i + 3, a] 
            out += line_form.format(atom, x, y, z, dx, dy, dz)
        out += '\n'

    with open('frequencies.xyz', 'w') as f:                     # writing output to file frequencies.xyz
        f.write(out)

# will run the following if calling directly (for testing purposes)
if __name__=='__main__':

    mol = Molecule.from_file('../extra-files/molecule.xyz')   #building molecule from molecule.xyz
    get_freqs(mol,'../extra-files/hessian.dat')
