#!/Users/avery/anaconda/bin/python

#####  Unit conversions
bohr2m = 5.2917721e-11
amu2kg = 1.6605389e-27
hartree2J = 4.3597443e-18
c = 29979245800.0 

import numpy as np
import sys
sys.path.insert(0, '../../0/aewiens')
from molecule import Molecule

class Frequencies:

    def __init__(self,mol,hessian_file):

        self.mol = mol
        self.hessian = hessian_file
        self.N = mol.__len__()

        m = []
        for i in range(self.N):
            m += [1/(mol.M[i])**0.5]*3
        self.MM = np.diag(m)
        self.m = m

    def get_MWhessian(self):

        h  = open(self.hessian,"r")
        H0 = np.matrix( [i.split() for i in h.readlines()],float )
        mwH = np.dot(self.MM, np.dot(H0,self.MM) )
        return mwH

    def get_frequencies(self):

        self.e, self.l = np.linalg.eigh( self.get_MWhessian() )
        self.Q = np.matrix(self.MM)*np.matrix(self.l)

        freq = []
        conv =  np.sqrt(hartree2J/(amu2kg*bohr2m*bohr2m))/(c*2*np.pi)  # dimensional analysis
        for i in self.e:
            if i <0:
                freq.append((-i)**0.5*conv)
            else:
                freq.append(i**0.5*conv)

        return freq

    def visualize_frequencies(self, output):

        mol = self.mol
        freq = self.get_frequencies()

        t = open(output,"w")
        for i in range(3*self.N):
            t.write("%d\n%s cm^{-1}\n" % (self.N, str(freq[i])))
            for j in range(self.N):
                atom = mol.atoms[j]
                x,y,z = mol.geom[j,0], mol.geom[j,1], mol.geom[j,2]
                dx,dy,dz = self.Q[3*j,i], self.Q[3*j+1,i], self.Q[3*j+2,i]
                #t.write("{:s}".format(atom) + "{0:20.12f}{1:20.12f}{1:20.12f}".format(x,y,z) )
                t.write("%s%15.7f%15.7f%15.7f%15.7f%15.7f%15.7f\n" % (atom, x, y, z,dx,dy,dz))
            t.write("\n")

        return None
        

f = open("../../extra-files/molecule.xyz", "r")
f = f.read()
mol = Molecule(f,"Bohr")

test = Frequencies(mol,"../../extra-files/hessian.dat")

test.visualize_frequencies("modes.xyz") 
