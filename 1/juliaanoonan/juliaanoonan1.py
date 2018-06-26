#Okay, let's do this

import sys
import math
import numpy as np
sys.path.insert(0, '../../extra-files')
import masses

def frequencies(mol,Hess):

#3 build the mass-weighted Hessian (do simple example to convince yourself)
    mas = np.array(np.repeat(mol.masses, 3) ** -0.5)
    mas = np.diag(mas)
    MwH = mas @ Hess @ mas

#4 compute the eigenvalues and eigenvectors of the mass-weighted Hessian

    k, q = np.linalg.eigh(MwH)

#5 un-mass-weight eigenvectors to get normal coords (makes more sense on paper) remember the eigenvalues same

    Q = mas @ q

#6 determine the spatial frequencies in wavenumbers from force constants (eigenvalues)

    htoJ = 4.359744e-18
    a0tom = 5.29177e-11
    utokg = 1.66054e-27
    c = 29979245800 #m/s
    tonewk = (htoJ)/((a0tom**2)*utokg)
    const = 1/(2*np.pi*c) #constants from finding frequencies and converting to wavenumbers

    knew = tonewk*k#gives force constant in right units
    v = []
    for val in knew:
        v = v + [const * np.sqrt(val.astype(complex))]

   # v = const*np.sqrt(knew.astype(complex))
    print(v)
#7 write normal modes to .xyz file, with frequencies neatly in comment lines

    with open(mol.xyz,'r') as f:
        lines = [line.rstrip() for line in f]
    out = ""
    for i in range(len(v)):
        out += lines[0]
        out += '\n'
        if np.isreal(v[i]):
            out += '{: >7.3f} cm^-1'.format(np.absolute(v[i]))
        else:
            out += '{: >7.3f}i cm^-1'.format(np.absolute(v[i]))
        out += '\n'
        for j in range(mol.natom): #0, 1, 2
            coord = lines[2+j].split()
            x = coord[1]
            y = coord[2]
            z = coord[3]
           # x, y, z = [coord][1], coord[2], coord[3]
            dx, dy, dz = Q[j:j+7:3,i]   #getting these from the eigenvectors in normal coords
            line = '{:<3}{:>15}{:>15}{:>15}{:>25.20f}{:>25.20f}{:>25.20f}'.format(mol.labels[j],x,y,z,dx,dy,dz)
            out += line+'\n'
        out += '\n'

    with open('frequencies.xyz', 'w+') as f:                     # write the xyz
        f.write(out)


#to test

if __name__ == '__main__':
    sys.path.insert(0, '../../0/juliaanoonan')
    from juliaanoonan0 import Molecule

#1 read in and build the Hessian

    Hesslinelist = open('../extra-files/hessian.dat','r').readlines()
    Hess = []
    for line in Hesslinelist:
        line = [float(i) for i in line.split()]
        Hess.append(line)
    Hess = np.array(Hess)

#2 build a Molecule object

    mol = Molecule('../../extra-files/molecule.xyz')

#call this bad boy in to get some displacements in that sweet, sweet xyz format
    frequencies(mol,Hess)
