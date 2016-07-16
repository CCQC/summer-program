import numpy as np
import sys
import os

sys.path.insert(0, '../../0/oliviabern')
from molecule import Molecule
sys.path.insert(0, '../../extra-files')
from masses import mass

#build molecule object
geom_string = open('../../extra-files/molecule.xyz').read()
mol = Molecule(geom_string)

atoms , xyz = mol.read(geom_string)
#xyz[0][0] += .005
m = str(mol)
m = m.splitlines()
N = int(m[0])


def inputform(atoms, xyz):
    line_form = '{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n'
    out = '{:d} {:d}\n'.format(0 , 1)
    for i in range(N):
        out += line_form.format(atoms[i] , *xyz[i]) 
    return out


#build input file template
template = open('../../extra-files/template.dat').read()
temp = 'memory 256 mb\n' + template
temp = temp.replace('{','',1)
temp = temp.replace('}','',2)

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):
    atoms , xyz = mol.read(geom_string)
    out = inputform(atoms, xyz)
    p = '/Users/oliviambernstein/Desktop/ccqc/'
    ph = p + 'hessian/'
    if not os.path.exists(ph):
        os.mkdir(ph)
    if not os.path.exists(ph + 'd/'):
        os.mkdir(ph + 'd/')
        t = temp.replace('{:s',out)
        file = open(ph + 'd/input.dat','w')
        file.write(t)
        file.close()
        """
    for i in range(N):
        for k in range(N):
            disp = xyz
            disp[k][i] += disp_size
            out = inputform(atoms, disp)
            name = ph + 'd+_' + str(k+1) + str(i+1) + '/'
            if not os.path.exists(name):
                os.mkdir(name)
                t = temp.replace('{:s' , out)
                file = open(name + 'input.dat', 'w')
                file.write(t)
                file.close()
            atoms , xyz = mol.read(geom_string)
    for i in range(N):
        for k in range(N):
            disp = xyz
            disp[k][i] -= disp_size
            out = inputform(atoms, disp)
            name = ph + 'd-_' + str(k+1) + str(i+1) + '/'
            if not os.path.exists(name):
                os.mkdir(name)
                t = temp.replace('{:s' , out)
                file = open(name + 'input.dat', 'w')
                file.write(t)
                file.close()
            atoms , xyz = mol.read(geom_string)
            """
    #positive displacements
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    if (i==k) and (j == l):
                        atoms , xyz = mol.read(geom_string)
                        disp = xyz
                        disp[i][j] += disp_size
                        out = inputform(atoms, disp)
                        name = ph + 'd+_' + str(i+1) + str(j+1) + '/'
                        if not os.path.exists(name):
                            os.mkdir(name)
                            t = temp.replace('{:s' , out)
                            file = open(name + 'input.dat', 'w')
                            file.write(t)
                            file.close()
                    else:
                        atoms , xyz = mol.read(geom_string)
                        disp = xyz
                        disp[i][j] += disp_size
                        disp[k][l] += disp_size
                        out = inputform(atoms,disp)
                        name = ph + 'd+_' + str(i+1) + str(j+1) + '_' + str(k+1) + str(l+1) + '/'
                        if not os.path.exists(name):
                            os.mkdir(name)
                            t = temp.replace('{:s' , out)
                            file = open(name + 'input.dat', 'w')
                            file.write(t)
                            file.close()
    #negative displacements
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    if (i==k) and (j == l):
                        atoms , xyz = mol.read(geom_string)
                        disp = xyz
                        disp[i][j] -= disp_size
                        out = inputform(atoms, disp)
                        name = ph + 'd-_' + str(i+1) + str(j+1) + '/'
                        if not os.path.exists(name):
                            os.mkdir(name)
                            t = temp.replace('{:s' , out)
                            file = open(name + 'input.dat', 'w')
                            file.write(t)
                            file.close()
                    else:
                        atoms , xyz = mol.read(geom_string)
                        disp = xyz
                        disp[i][j] -= disp_size
                        disp[k][l] -= disp_size
                        out = inputform(atoms,disp)
                        name = ph + 'd-_' + str(i+1) + str(j+1) + '_' + str(k+1) + str(l+1) + '/'
                        if not os.path.exists(name):
                            os.mkdir(name)
                            t = temp.replace('{:s' , out)
                            file = open(name + 'input.dat', 'w')
                            file.write(t)
                            file.close()
                        """
    for i in range(N):
        for j in range(N):
            disp = xyz
            disp[i][j] -= disp_size
            for k in range(N):
                for l in range(N):
                    if (i==k) and (j == l):
                        print i    
                    else:
                        disp[k][l] -= disp_size
                        out = inputform(atoms,disp)
                        name = ph + 'd-_' + str(i+1) + str(j+1) + '_' + str(k+1) + str(l+1) + '/'
                        if not os.path.exists(name):
                            os.mkdir(name)
                            t = temp.replace('{:s' , out)
                            file = open(name + 'input.dat', 'w')
                            file.write(t)
                            file.close()
                        atoms , xyz = mol.read(geom_string)
                        """
                    



generate_inputs(mol, template)
