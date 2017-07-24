# coding: UTF-8

import math
import numpy as np
import sys
sys.path.insert(0, '../../0/misiewicz')
from molecule import Molecule

def frequencies(mol, hessian):
    '''Takes a molecule object and a Hessian matrix (Hermitian n-by-n matrix)
    and writes the vibrational frequencies and modes to normal_modes.xyz.'''
    
    if not isinstance(mol, Molecule):
        raise Exception("Please pass a molecule to frequencies.")
    if not np.array_equal(hessian, np.matrix.getH(hessian)):
        raise Exception("Please pass a Hermitian matrix to frequencies.")
    if hessian.shape != (3*mol.natom, 3*mol.natom):
        raise Exception("Dimensions of the molecule and the Hessian disagree.")
    mol.to_angstrom()

    # A 3n array used for massweighting, per eqn. 2.
    # A 3n-by-3n matrix could be made via broadcasting, but this feels cleaner.
    massweighter = np.repeat(np.asarray(mol.masses)**(-0.5), 3)
    # Multiply the columns, (the B terms) and then the rows (the A terms).
    # Normally, we would need to un-transpose back, but hessian is Hermitian.
    hessian = np.transpose(hessian * massweighter) * massweighter
    evalues, evectors = np.linalg.eigh(hessian)
    # Eigenvectors are the columns of evectors.
    # q_ba is the bth element of the ath eigenvector
    # Undo the massweight.
    evectors = np.transpose(evectors.T * massweighter)
    # Convert units to S.I.
    evalues *= (4.3597443e-18/5.2917721e-11**2/1.6605389e-27)
    frequencies = 1/(np.pi*2)*np.sqrt(evalues.astype(complex))/29979245800

    file_contents = ""
    format_str = "{:3}" + " {:15.10f}" * 6 + "\n"

    for i, frequency in enumerate(frequencies):
        if np.isreal(frequency):
            freq_prnt = str(np.real(frequency)) + " cm^-1"
        else:
            freq_prnt = str(frequency) + " Complex!"
        file_contents += (str(mol.natom) + "\n" + str(freq_prnt) + "\n")
        # The eigenvectors are columns, so pull from index 2.
        mode = evectors[:, i].reshape(mol.natom, 3)
        for label, geom, atom in zip(mol.labels, mol.geom.tolist(), mode):
            file_contents += format_str.format(label, *geom, *atom)
        file_contents += "\n\n"

    with open("normal_modes.xyz", "w") as f:
        f.write(file_contents)

if __name__ == "__main__":
    hessian = np.loadtxt("../../extra-files/hessian.dat")
    water = Molecule("../../extra-files/molecule.xyz")
    frequencies(water, hessian)