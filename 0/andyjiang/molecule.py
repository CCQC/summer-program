import sys
import numpy as np
sys.path.append('../../extra-files')
from masses import *

BOHR_TO_ANGSTROM = 0.529177
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM

class Molecule:
    def __init__(self, contents : str):
        lines = contents.split("\n")
        self.natom = int(lines[0])
        self.units = lines[1]
        self.labels = [None] * self.natom
        self.charges = [None] * self.natom
        self.geoms = np.zeros((self.natom, 3))

        for idx in range(self.natom):
            pieces = lines[idx+2].split()
            self.labels[idx], self.geoms[idx][0], self.geoms[idx][1], self.geoms[idx][2] = \
                pieces[0].upper(), float(pieces[1]), float(pieces[2]), float(pieces[3])
            self.charges[idx] = get_charge(pieces[0])

    def __len__(self):
        return self.natom

    def __str__(self):
        mol_str = f"Molecule info:\nNumber of Atoms: {self.natom}\nElements:  "
        for label in self.labels:
            mol_str += f"{label}  "
        return mol_str
    
    def __repr__(self):
        mol_str = f'{str(self.natom)}\n{self.units}'
        for idx in range(self.natom):
            mol_str += f'\n{self.labels[idx]}    {self.geoms[idx][0]:.8f}    {self.geoms[idx][1]:.8f}    {self.geoms[idx][2]:.8f}'
        return mol_str

    def __iter__(self):
        self.idx = 0
        return self
    
    def __next__(self):
        if (self.idx <= self.natom - 1):
            value = (self.labels[self.idx], self.geoms[self.idx][0], self.geoms[self.idx][1], self.geoms[self.idx][2])
            self.idx += 1
            return value
        else:
            raise StopIteration

    def __add__(self, mol):

        natom = self.natom + mol.natom
        units = self.units

        if (units == "Angstrom"):
            mol.to_angstrom()
        else:
            mol.to_bohr()

        mol_str = f'{str(natom)}\n{units}'

        for idx in range(self.natom):
            mol_str += f'\n{self.labels[idx]}    {self.geoms[idx][0]:.8f}    {self.geoms[idx][1]:.8f}    {self.geoms[idx][2]:.8f}'
        
        for idx in range(mol.natom):
            mol_str += f'\n{mol.labels[idx]}    {mol.geoms[idx][0]:.8f}    {mol.geoms[idx][1]:.8f}    {mol.geoms[idx][2]:.8f}'

        return mol_str

    def to_bohr(self):
        if self.units == "Angstrom":
            self.units = "Bohr"
            self.geoms *= ANGSTROM_TO_BOHR
    
    def to_angstrom(self):
        if self.units == "Bohr":
            self.units = "Angstrom"
            self.geoms *= BOHR_TO_ANGSTROM

    def xyz_string(self):
        return repr(self)

    def copy(self):
        return Molecule(self.xyz_string())

if __name__ == '__main__':
    with open('../../extra-files/molecule.xyz', 'r') as mol_file:
        contents = mol_file.read()

    mol = Molecule(contents)
    print(mol.natom)
    print(mol.units)
    print(mol.labels)
    print(mol.charges)
    print(mol.geoms)

    mol.to_bohr()
    print(mol.units)
    print(mol.geoms)
    
    mol.to_angstrom()
    print(mol.units)
    print(mol.geoms)

    mol2 = mol.copy()
    print(mol2.natom)
    print(mol2.units)
    print(mol2.labels)
    print(mol2.charges)
    print(mol2.geoms)

    for element in mol2:
        print(element)

    print(repr(mol2))
    print(len(mol2))
    print(str(mol2))

    mol.to_bohr()

    print(mol + mol2)