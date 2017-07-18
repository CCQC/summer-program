import frequencies
import geom

hessian = frequencies.readhessian('hessian.dat')
molecule = geom.Molecule('molecule.xyz')

vecs, vals = frequencies.frequencies(molecule,hessian)
for idx, val in enumerate(vals):
     print(val, vecs[idx])
