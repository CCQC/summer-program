import sys
sys.path.append('../../0/mdav2')
sys.path.append('../../extra-files')
import frequencies
import molecule

#TODO: clean this up into a better test
#TODO: learn how to test
hessian = frequencies.readhessian('../../extra-files/hessian.dat')
molecule = molecule.Molecule('../../extra-files/molecule.xyz')

vecs, vals = frequencies.frequencies(molecule,hessian)

frequencies.printfrequencies(molecule,vals,vecs)
