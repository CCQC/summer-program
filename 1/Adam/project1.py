import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/Adam')

from molecule import Molecule

mol = Molecule(open("../extra-files/molecule.xyz").read())

mol = mol.bohr_to_ang()
#read Hessian Values
hessian_values = open('../extra-files/hessian.dat').readlines()
print (hessian_values)
#collect hessain values into matrix
for row in hessian_values:
    row.split()






