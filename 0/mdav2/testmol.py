import molecule
x = molecule.Molecule('../../extra-files/molecule.xyz')
print(x.xyz_string())
y = x.copy()
print(y.xyz_string())
#Write unittests etc
