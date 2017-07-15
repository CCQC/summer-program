import frequencies
import geom

hessian = frequencies.readhessian('hessian.dat')
molecule = geom.Molecule('molecule.xyz')

frequencies.frequencies(molecule,hessian)
