#! /usr/bin/python4

import numpy as np
import math as m
import cmath as cm
import sys
sys.path.insert(0, '../extra-files')
sys.path.insert(0, '../../0/aroeira')
from molecule import Molecule
from masses import get_mass

c = 2.998e10
E = 4.359744e-18
a0 = 5.291772e-11
u = 1.66054e-27

def mult(a, b):
    b = np.transpose(b)
    l = len(a)
    r = np.zeros((l, l))
    for i, roll in enumerate(a):
        for j, collum in enumerate(b):
            r[i][j] = np.dot(roll, collum)
    return r

class Hessian:
    def __init__(self, molecule, matrix):
        self.mol = molecule
        self.read(matrix) 
        self.frequencies()

    def __str__(self):
        out = ''
        d = np.transpose(self.disp).tolist()
        for w, wnum in enumerate(self.wavenumbers):
            out += str(len(self.mol)) + '\n'
            if wnum.imag > 0:
                out += '{:<0.2f}i cm^-1\n'.format(wnum.imag)
            else:
                out += '{:<5.2f} cm^-1\n'.format(wnum.real)
            for a, atom in enumerate(self.mol.atoms):
                out += '{:<3.8s}'.format(atom)
                for c in self.mol.geom[a]:
                    out += '{:>20.13f}'.format(c)
                for x in d[w][:3]:
                    out += '{:>20.13f}'.format(x)
                d[w] = d[w][3:]
                out += '\n'
            out += '\n'    
        return out

    def read(self, matrix):
        with open(matrix) as lines:
            self.hessian = []
            for line in lines:
                hold = []
                for term in line.split():
                    hold.append(float(term))
                self.hessian.append(hold)
            self.hessian = np.array(self.hessian)

    def masses(self):
        masses = []
        for atom in self.mol.atoms:
            masses.append(1/np.sqrt(get_mass(atom)))
            masses.append(1/np.sqrt(get_mass(atom)))
            masses.append(1/np.sqrt(get_mass(atom)))
        M = np.diag(masses)
        return M

    def weighted(self):
        mw_hessian = mult(mult(self.masses(), self.hessian), self.masses())
        return mw_hessian

    def frequencies(self):
        freq, coord =  np.linalg.eigh(self.weighted())
        self.disp = mult(self.masses(), coord)
        self.wavenumbers = []
        for r in freq:
            self.wavenumbers.append(cm.sqrt(r*E/(a0*a0*u))/(2*np.pi*c))


    def output(self):
        f = open('output.xyz', 'w')
        f.write(str(self))
        f.close()    

water = Molecule('../extra-files/molecule.xyz','Bohr')
hessian = Hessian(water, '../extra-files/hessian.dat')
hessian.output()
