#!/usr/bin/env python

import numpy as np
import os
import re
os.sys.path.insert(0,"../../0/aewiens")

from molecule import Molecule


class Hessian(object):
    """
    Compute a finite difference Hessian
    """

    def __init__(self,mol,template,disp_size=0.005):

        self.mol = mol
        self.template_str = open(template,"r").read()
        self.N = len(self.mol)
        self.h = disp_size

    def make_input(self,dirname,atoms,geom):
        """
        :param dirname: directory for writing a psi4 input file
        :param atoms: list of strings of atom labels
        :param geom: 2darray of xyz coordinates
        """
        os.mkdir("%s" % dirname)
        xyz = [('  ').join([atoms[i]] + [str(geom[i,j]) for j in range(3)]) for i in range(self.N)]
        open("{:s}/input.dat".format(dirname), 'w').write(self.template_str.format(('\n').join(xyz)))


    def run_input(self,dirname):
        """
        Run a psi4 single-point energy with specified name
        :param dirname: directory where we'll run the input file
        """
        os.chdir(dirname)
        os.system('psi4')
        os.chdir('..')
    
    
    def find_E(self,i,j,hi,hj):
        """
        Pull energy from output.dat and throw an error if not found
        :params i,j: indices of atoms 1,2
        :params hi,hj: displacements of atoms 1,2 (-1, 0, or 1, corresponds to -h, 0, or h)
        """
        dirname = "X%dX%d_%d%d" % (i,j,hi,hj)
        out_str = open("disps/%s/output.dat" % dirname, "r").read()
        match = re.findall("Total Energy\s=\s+-\d+.\d+",out_str)
        if match == []:
            out = "Cannot find energy!"
        else:
            out = float(match[0].split()[-1])
        return out


    def run_disps(self):

        h, N = self.h, self.N
        os.mkdir("disps")
        os.chdir("disps")

        self.make_input("X0X0_00",mol.atoms,mol.geom)
        self.run_input("X0X0_00")


        ####   Run single displacements   ####

        for i in range(3*N):
           forward = "X%dX0_10" % i
           reverse = "X%dX0_-10" % i
           geom_copy = mol.copy().geom
           geom_copy[i/3,i%3] +=h

           self.make_input(forward,mol.atoms,geom_copy)
           self.run_input(forward)

           geom_copy[i/3,i%3] -= 2*h
           self.make_input(reverse,mol.atoms,geom_copy)

           self.run_input(reverse)

        ####   Run double displacements    ######

        for i in range(3*N):
            for j in range(i):
                forward = "X%dX%d_11" % (i,j)
                reverse = "X%dX%d_-1-1" % (i,j)
                geom_copy2 = mol.copy().geom

                geom_copy2[i/3,i%3] += h
                geom_copy2[j/3,j%3] += h

                self.make_input(forward,mol.atoms,geom_copy2)

                geom_copy2[i/3,i%3] -= 2*h
                geom_copy2[j/3,j%3] -= 2*h

                self.make_input(reverse,mol.atoms,geom_copy2)

                self.run_input(forward)
                self.run_input(reverse)

        os.chdir("..")


    def make_Hessian(self):

        self.run_disps()
        
        h, N = self.h, self.N
        E0 = self.find_E(0,0,0,0)
        self.H = np.zeros((3*self.N, 3*self.N))

        for i in range(3*N):
            for i in range(3*N):
                self.H[i,i]= (self.find_E(i,0,1,0)+self.find_E(i,0,-1,0)-2*E0)/(h**2)
                for j in range(0,i):
                    self.H[i,j] = (self.find_E(i,j,1,1)+self.find_E(i,j,-1,-1)-self.find_E(i,0,1,0)-self.find_E(j,0,1,0)-self.find_E(j,0,-1,0)-self.find_E(i,0,-1,0)+2*E0)
                    self.H[i,j] /= 2*h**2
                    self.H[j,i] = self.H[i,j]


    def write_Hessian(self):
        """
        write Hessian matrix to hessian.dat file
        """
        self.make_Hessian()
        np.savetxt("hessian.dat",self.H,"%15.7f"," ","\n")


if __name__ == "__main__":

    mol = Molecule(open("../../extra-files/molecule.xyz","r").read() )
    mol.bohr()
    hessian = Hessian(mol,"template.dat")
    hessian.write_Hessian()
