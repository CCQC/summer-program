import sys
sys.path.insert(0, '../../0/bz3wp/')
sys.path.insert(0, '../extra-files')
import masses
from molecule import Molecule
import numpy as np

#building molecule from molecule.xyz
mol = Molecule.from_file('../extra-files/molecule.xyz')
mol.to_angstrom()

#reading in hessian matrix as numpy array 
hessian = np.loadtxt('../extra-files/hessian.dat')

def frequencies(mol, hessian):
    mass = []
    for atom in mol.labels:
        mass += [masses.get_mass(atom)**(-0.5)]*3
    M = np.diag(mass) 
    
    #building mass weighted hessian 
    H = M @ hessian @ M
    
    #diagonalizing weighted hessian 
    k, q = np.linalg.eigh(H)

    #un-mass-weight eigenvectors
    Q = M @ q

    #determining spatial frequencies in cm^-1
    hartreetoJ = 4.3597443e-18
    amutokg = 1.6605389e-27
    bohrtom = 5.2917721e-11
    c = 2.99792458e10
    ka = k * hartreetoJ * (1/amutokg) * (1/bohrtom)**2
    va_conv = (1/(2 * np.pi)) * (1/c)
    va = ka * va_conv
      
    #formatting for output  
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
            dx, dy, dz = q[3*i: 3*i + 3, a] 
            out += line_form.format(atom, x, y, z, dx, dy, dz)
        out += '\n'

    with open('frequencies.xyz', 'w') as f:
        f.write(out)
 
frequencies(mol, hessian)
