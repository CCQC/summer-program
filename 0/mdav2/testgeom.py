import geom
x = geom.Molecule('molecule.xyz')
print(x.geom)
x.to_bohr()
print(x.geom)
x.to_angstrom()
print(x.geom)
