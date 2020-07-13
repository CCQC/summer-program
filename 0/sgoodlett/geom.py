import masses
import numpy as np
from copy import deepcopy

class molecule():

    def __init__(self, filename='molecule.xyz', units='Angstrom',charge=0,mult=1):
        self.units = units
        self.labels = []
        self.charges = []
        self.masses = []
        self.geom_list = []
        self.index = 0
        self.charge = charge
        self.mult = mult
        with open(filename,'r') as fn:
            self.natom = int(fn.readline())
            geom_lines = fn.readlines()[1:]
            for line in geom_lines:
                line_list = line.split()
                self.labels.append(line_list[0].upper())
                self.charges.append(masses.get_charge(line_list[0]))
                self.masses.append(masses.get_mass(line_list[0]))
                self.geom_list.append([line_list[1],line_list[2],line_list[3]])
        self.geom = np.array(self.geom_list).astype(float)

    def to_bohr(self):
        if self.units != 'Bohr':
            self.geom = self.geom * 1.88973
            self.units = 'Bohr'
        else:
            print('Already in Bohr')

    def to_angstrom(self):
        if self.units != 'Angstrom':
            self.geom = self.geom * 0.529177
            self.units = 'Angstrom'
        else:
            print('Already in Angstroms')

    def xyz_string(self,option=0):
        if option == 0:
            new_string = str(self.natom) + '\n'
        else:
            new_string = '{} {}\n'.format(self.charge,self.mult)
        for i in range(len(self.labels)):
            new_line = '{} {} {} {}'.format(self.labels[i],self.geom[i][0],self.geom[i][1],self.geom[i][2])
            new_string = new_string + new_line + '\n'
        return new_string

    def copy(self):
        return deepcopy(self)

    def __len__(self):
        return self.natom

    def __str__(self):
        new_string = str(self.natom) + self.units + '\n'
        for i in range(len(self.labels)):
            new_line = '{} {} {} {} {} {}'.format(self.labels[i],self.charges[i],self.masses[i],self.geom[i][0],self.geom[i][1],self.geom[i][2])
            new_string = new_string + new_line + '\n'
        return new_string

    def __repr__(self):
        return self.xyz_string()

    def __iter__(self):
        return self

    def __next__(self):
        if self.index >= self.natom:
            raise StopIteration
        self.index += 1
        return (self.labels[self.index-1],self.geom[self.index-1])

    def __add__(self, other):
        mol = self.copy()
        if self.units != 'Angstrom':
            self.to_angstrom()
        if other.units != 'Angstrom':
            other.to_angstrom()
        mol.geom = np.concatenate((self.geom,other.geom),axis=0)
        mol.labels = self.labels + other.labels
        mol.charges = self.charges + other.charges
        mol.masses = self.masses + other.masses
        mol.natom = self.natom + other.natom
        return mol

if __name__ == "__main__":
    m = molecule()
    print(m.xyz_string())
    m.to_bohr()
    print(m.xyz_string())
    #Test copy, len, str, repr, iter, add
    n = m.copy()
    print(len(m))
    print(str(m))
    for i in m:
        print(i)
    a = m + n
    print(a.xyz_string())
