import sys
sys.path.append('../../0/zachglick')

import math
from numpy import matrix
from numpy.linalg import eigh
from molecule import Molecule

# Some physical constants. Maybe I should make a seperate class for these?
joule_per_hartree = 4.35974e-18
bohr_per_meter = 5.29177e-11
kilogram_per_amu = 1.6605389e-27
centimeter_per_second = 2.99792458e10

def format_frequency(eigval) :
    """ 
    eigval -- an eigenvalue of some hessian in units of Hartree / (bohr * bohr * kilogram)
    
    Returns (str) the frequency in cm^-1 represented by the eigenvalue (formatted nicely)
    Note that the frequency can be imaginary
    """
    
    # To avoid working with complex numbers
    imag_freq = (eigval < 0.0)
    eigval = abs(eigval)

    # Convert to rad^2 / sec^2
    eigval = eigval * joule_per_hartree / ( (bohr_per_meter ** 2) * kilogram_per_amu )
    
    # Then to Hz
    eigval = math.sqrt(eigval) / (2 * math.pi * centimeter_per_second)
    
    if imag_freq :
        return "{:9.3f}i cm^-1".format(eigval)
    else :
        return "{:9.3f}  cm^-1".format(eigval)

def output_frequencies(mol, hess) :
    """
    mol -- an object of type Molecule
    mat -- a Hessian matrix for the Molecule
    
    Calculated (and writes to frequencies.out) the associated vibrational modes 
    """
    
    # Iterating through the hessian
    for row_ind in range(mol.natom*3):
        for col_ind in range(mol.natom*3):
            
            # Get the atom index from the location in the hessian
            row_atom = row_ind // 3
            col_atom = col_ind // 3
            
            # Mass weighting
            if row_atom == col_atom :
                hess[row_ind, col_ind] /= mol.masses[row_atom] # sqrt(x ** 2) = x
            else :
                hess[row_ind, col_ind] /= math.sqrt(mol.masses[row_atom]*mol.masses[col_atom])

    # Solve with eigh() instead of eig() because hessian is symmetric
    eigvals, eigvecs = eigh(hess)

    # Unweight the eigenvectors
    for atom_dim_ind in range(len(eigvecs)):
        atom_mass = mol.masses[atom_dim_ind // 3]
        eigvecs[atom_dim_ind] /= math.sqrt(atom_mass)
    # Write results to file for display in jmol
    with open('frequency.out', 'w') as f:
        for freq_ind in range(3*mol.natom) :
            f.write( "{}\n".format(mol.natom) )
            f.write( format_frequency(eigvals[freq_ind]) + "\n" )
            for atom_ind in range(mol.natom) :
                
                # Get and format the initial coordinates from the Molecule coordinates
                coords = "".join([ "{:>16.12f}".format(coord) for coord in mol.geom[atom_ind]])
                
                #Get and format the displacement from the eigenvector
                displs = eigvecs[3*atom_ind:3*(atom_ind+1),freq_ind]
                displs = "".join([ "{:>16.12f}".format(disp[0,0]) for disp in displs])
                
                f.write( "%s %s %s\n" % (mol.labels[atom_ind], coords, displs) )
            f.write('\n')

def solve(mol_path, hess_path):

    with open(mol_path, 'r') as f:
        molecule = Molecule(f.read())
        molecule.to_angstrom()

    with open(hess_path, 'r') as f:
        str = (f.read()).replace("\n",";")
        while str[-1] == ';' :
            str = str[:-1]
        mat = matrix(str)

    output_frequencies(molecule, mat)

if __name__ == '__main__':
    
    solve('../../extra-files/molecule.xyz', '../../extra-files/hessian.dat')

