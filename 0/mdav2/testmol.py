import molecule
x = molecule.Molecule('../../extra-files/molecule.xyz')
print('The original file')
print(x.xyz_string())

print('Copy of the molecule object')
y = x.copy()
print(y.xyz_string())
#Write unittests etc
