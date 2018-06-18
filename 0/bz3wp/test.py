# Tests molecule.py by calling static method on molecule.xyz

from molecule import Molecule

mol=Molecule.from_file('../../extra-files/molecule.xyz')
print(mol)
