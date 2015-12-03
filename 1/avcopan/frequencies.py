import sys
sys.path.insert(0, "../../0/avcopan/")

from molecule import Molecule
import numpy as np

# 1. read in Hessian matrix
H     = np.loadtxt("../../extra-files/hessian.dat")

# 2. read in molecule object
mol   = Molecule(open("../../extra-files/molecule.xyz").read())

# 3. build mass-weighted Hessian matrix
natom = mol.natom
m     = np.repeat(mol.masses, 3) ** -0.5
M     = np.broadcast_to(m, (3*natom, 3*natom)).T
tH    = M.T * H * M

# 4. compute eigenvalues and eigenvectors of the mass-weighted Hessian matrix
k, tQ  = np.linalg.eigh(tH)

# 5. un-mass-weight the normal coordinates
Q = M * tQ

# 6. get wavenumbers from a.u. force constants
hartree2J = 4.3597443e-18
amu2kg = 1.6605389e-27
bohr2m = 5.2917721e-11
c = 29979245800.0 # speed of light in cm/s
v = np.sqrt(np.complex_(k)) * np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)

ret = ''
fmt = '{:2s}' + '{: >15.10f}'*6 + '\n'
for A in range(3*natom):
  if np.isreal(v[A]):
    ret += '{:d}\n{: >7.2f}  cm^-1\n'.format(natom, np.absolute(v[A]))
  else:
    ret += '{:d}\n{: >7.2f}i cm^-1\n'.format(natom, np.absolute(v[A]))
  for B in range(natom):
    label      = mol.labels[B]
    x ,  y,  z = mol.geom[B]
    dx, dy, dz = Q[3*B:3*B+3, A]
    ret += fmt.format(label, x, y, z, dx, dy, dz)
  ret += '\n'

print ret

# 7. save the normal modes to a file
outfile = open('modes.xyz', 'w+')
outfile.write(ret)
outfile.close()
