#! /usr/bin/python3

import numpy as np
import math as m
import cmath as cm
from molecule import Molecule
from masses import get_mass

c = 2.998e10
E = 4.359744e-18
a0 = 5.291772e-11
u = 1.66054e-27

def get_hessian(h_file):
    with open(h_file) as lines:
        hessian = []
        for line in lines:
            hold = []
            for term in line.split():
                hold.append(float(term))
            hessian.append(hold)
        hessian = np.array(hessian)
    return hessian

def mass_weighted(hessian, molecule):
    d = len(hessian)
    mw_hessian = []
    i_index = 0
    for i in range(d):
        hold = []
        j_index = 0
        if i%3 == 0 and i != 0:
            i_index += 1
        atom_i = molecule.atoms[i_index]
        for j in range(d):
            if j%3 == 0 and j != 0:
                j_index += 1
            atom_j = molecule.atoms[j_index]
            w = m.sqrt(get_mass(atom_i)*get_mass(atom_j))
            hold.append(hessian[i][j] / w)
        mw_hessian.append(hold)
    mw_hessian = np.array(mw_hessian)
    return mw_hessian    

def uneig(evec, molecule):
    d = len(evec)
    out = []
    index = 0
    for i in range(d):
        hold = []
        if i%3 == 0 and i != 0:
            index += 1
        atom = molecule.atoms[index]
        w = m.sqrt(get_mass(atom))
        for j in range(d):
            hold.append(evec[i][j] / w)
        out.append(hold)
    out = np.array(out)
    return out    

def output(molecule, wnum, disp):
    out = ''
    disp = np.transpose(disp)
    for w_index in range(len(wnum)):
        d_index = 0
        out += str(len(molecule)) + '\n'
        if wnum[w_index].imag > 0:
            out += '{:<0.2f}i cm^-1\n'.format(wnum[w_index].imag)
        else:
            out += '{:<5.2f} cm^-1\n'.format(wnum[w_index].real)
        for a_index in range(len(molecule.atoms)):
            out += '{:<3.8s}'.format(molecule.atoms[a_index])
            for c_index in range(3):
                out += '{:>20.13f}'.format(molecule.geom[a_index][c_index])
            for x in range(3):
                out += '{:>20.13f}'.format(disp[w_index][d_index])
                d_index += 1
            out += '\n'
        out += '\n'
    return out

def frequencies(molecule, hessian):
    mw_hessian = mass_weighted(hessian, water)
    freq, coord =  np.linalg.eigh(mw_hessian)
    disp = uneig(coord, water)
    wavenumbers = []
    for r in freq:
        wavenumbers.append(cm.sqrt(r*E/(a0*a0*u))/(2*np.pi*c))
    f = open('output.xyz', 'w')
    f.write(output(water, wavenumbers, disp))
    f.close()    

water = Molecule('water.xyz',)
water.bohr()
hessian = get_hessian('hessian.dat')
frequencies(water, hessian)
