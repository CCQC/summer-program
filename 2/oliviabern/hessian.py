#!/usr/bin/python

import numpy as np
import sys
import os
import subprocess

sys.path.insert(0, '../../0/oliviabern')
from molecule import Molecule
sys.path.insert(0, '../../extra-files')
from masses import mass

#build molecule object
geom_string = open('../../extra-files/molecule.xyz').read()
mol = Molecule(geom_string)
atoms, xyz = mol.atoms, mol.xyz
mol.converttoB()
#atoms, xyz = mol.atoms, mol.xyz

#atoms , xyz = mol.read(geom_string)
m = str(mol)
m = m.splitlines()
N = int(m[0])


def inputform(atoms, xyz):
    line_form = '{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n'
    out = '{:10s}\n'.format('units Bohr')
    for i in range(N):
        out += line_form.format(atoms[i] , *xyz[i]) 
    return out


#build input file template
template = open('../../extra-files/template.dat').read()
temp = 'memory 256 mb\n' + template
temp = temp.replace('{','',1)
temp = temp.replace('}','',2)

geo = xyz.reshape((1,3*N))
geo = geo[0][:]

def generate_inputs(mol, template, disp_size = 0.005, directory = "DISPS"):
    out = inputform(atoms, xyz)
    ph = directory + '/'
    if not os.path.exists(ph):
        os.mkdir(ph)
    if not os.path.exists(ph + 'd/'):
        os.mkdir(ph + 'd/')
        t = temp.replace('{:s',out)
        file = open(ph + 'd/input.dat','w')
        file.write(t)
        file.close()
    #positive single displacements
    for i in range(3*N):
        disp = list(geo)
        disp[i] += disp_size
        disp = np.array(disp)
        disp = disp.reshape((3,N))
        out = inputform(atoms, disp)
        name = ph + 'd+' + str(i) + '/'
        if not os.path.exists(name):
            os.mkdir(name)
            t = temp.replace('{:s' , out)
            file = open(name + 'input.dat', 'w')
            file.write(t)
            file.close()
    #negative single displacements
    for i in range(3*N):
        disp = list(geo)
        disp[i] -= disp_size
        disp = np.array(disp)
        disp = disp.reshape((3,N))
        out = inputform(atoms, disp)
        name = ph + 'd-' + str(i) + '/'
        if not os.path.exists(name):
            os.mkdir(name)
            t = temp.replace('{:s' , out)
            file = open(name + 'input.dat', 'w')
            file.write(t)
            file.close()
    #positive double displacements
    for i in range(3*N):
        for j in range(i+1,3*N):
            disp = list(geo)
            disp[i] += disp_size
            disp[j] += disp_size 
            disp = np.array(disp)
            disp = disp.reshape((3,N))
            out = inputform(atoms, disp)
            name = ph + 'd+' + str(i) + str(j) + '/'
            if not os.path.exists(name):
                os.mkdir(name)
                t = temp.replace('{:s' , out)
                file = open(name + 'input.dat', 'w')
                file.write(t)
                file.close()
    #negative double displacements
    for i in range(3*N):
        for j in range(i+1,3*N):
            disp = list(geo)
            disp[i] -= disp_size
            disp[j] -= disp_size 
            disp = np.array(disp)
            disp = disp.reshape((3,N))
            out = inputform(atoms, disp)
            name = ph + 'd-' + str(i) + str(j) + '/'
            if not os.path.exists(name):
                os.mkdir(name)
                t = temp.replace('{:s' , out)
                file = open(name + 'input.dat', 'w')
                file.write(t)
                file.close()



def run_jobs(mol , command = 'psi4' , directory = 'DISPS'):
    ph ='/Users/oliviambernstein/git/summer-program/2/oliviabern/' + directory + '/'
    os.chdir(ph + 'd/')
    subprocess.call(command)
    for i in range(3*N): 
        name = 'd+' + str(i) + '/'
        os.chdir(ph + name)
        subprocess.call(command)
    for i in range(3*N):
        name = ph + 'd-' + str(i) + '/'
        os.chdir(name)
        subprocess.call(command)
    for i in range(3*N):
        name = ph + 'd+' + str(i) + '/'
        os.chdir(name)
        subprocess.call(command)
    for i in range(3*N):
        for j in range(i+1,3*N):
            name = ph + 'd+' + str(i) + str(j) + '/'
            os.chdir(name)
            subprocess.call(command)
    for i in range(3*N):
        for j in range(i+1,3*N):
            name = ph + 'd-' + str(i) + str(j) + '/'
            os.chdir(name)
            subprocess.call(command)

def grab_energy(mol, energy_prefix = '@DF-RHF Final Energy:'):
    f = open('output.dat').read()
    x = f.find(energy_prefix)
    x += len(energy_prefix)
    x += 3
    z = '' 
    for i in range(x,len(f)):
        if f[i] == '\n':
            break
        else:
            z += f[i]
    z = float(z)
    return z
    




def build_hessian(mol, energy_prefix = '@DF-RHF Final Energy:', disp_size = 0.005, directory = "DISPS"):
    ph ='/Users/oliviambernstein/git/summer-program/2/oliviabern/' + directory + '/'
    os.chdir(ph + 'd/')
    d = grab_energy(mol, energy_prefix)
    hess = np.zeros((3*N,3*N))
    for i in range(3*N):
        for k in range(3*N):
            if i ==k:
                os.chdir(ph + 'd+' + str(i))
                dp = grab_energy(mol)
                os.chdir(ph + 'd-' + str(i))
                dm = grab_energy(mol)
                h = (dp + dm - 2*d)/(disp_size**2)
                hess[i,k] += h
            else:
                j = sorted([i,k])
                os.chdir(ph + 'd+' + str(j[0]) + str(j[1]))
                dpp = grab_energy(mol)
                os.chdir(ph + 'd-' + str(j[0]) + str(j[1]))
                dmm = grab_energy(mol)
                os.chdir(ph + 'd+' + str(i))
                dpn = grab_energy(mol)
                os.chdir(ph + 'd-' + str(i))
                dmn= grab_energy(mol)
                os.chdir(ph + 'd+' + str(k))
                dnp = grab_energy(mol)
                os.chdir(ph + 'd-' + str(k))
                dnm = grab_energy(mol)
                h = (dpp + dmm - dpn -dmn - dnp - dnm + 2*d)/(2*disp_size**2)
                hess[i,k] += h
    hess = hess.tolist()

    out = ''
    f = ''
    f += '{: >12.7f} '*(3*N)
    f += ' \n'
    for i in range(3*N):
        out += f.format(*hess[i])

    os.chdir('/Users/oliviambernstein/git/summer-program/2/oliviabern')
    file = open('hessian.dat','w')
    file.write(out)
    file.close()



                
    

generate_inputs(mol, template)
#run_jobs(mol)
build_hessian(mol)
