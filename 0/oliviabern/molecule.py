import numpy as np

class Molecule:
    def __init__(self, geom_string, units = 'Angstrom'):
        self.read(geom_string)
        self.units = units

    def read(self , geom_string):
        lines = geom_string.splitlines()
        self.atoms = []
        xyz = []
        natoms = int(lines[0])
        units = lines[1]
        for line in lines[2:]:
            atom, x, y, z = line.split()
            self.atoms.append(atom)
            xyz.append([float(x),float(y), float(z)])
        self.xyz = np.array(xyz)
        return self.atoms , self.xyz

    def __len__(self):
        return len(self.atoms)

    def __str__(self):
        line_form = '{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n'
        out = '{:d}\n{:s}\n'.format(len(self) , self.units)
        for i in range(len(self)):
            out += line_form.format(self.atoms[i] , *self.xyz[i]) 
        return out

    def inputform(self):
        line_form = '{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n'
        out = '{:d} {:d}\n'.format(0 , 1)
        for i in range(len(self)):
            out += line_form.format(self.atoms[i] , *self.xyz[i]) 
        return out

    def converttoB(self):
        if self.units == 'Angstrom':
            self.xyz *= 1.889725989
            self.units = 'Bohr'


    def converttoA(self):
        if self.units == 'Bohr':
            self.xyz /= 1.889725989
            self.units = 'Angstrom'

if __name__ == "__main__":
    geom_str = open('../../extra-files/molecule.xyz').read()
    mol = Molecule(geom_str)
    mol.converttoB()
    mol.converttoA()
    print(mol)

