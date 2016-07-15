import numpy as np
from numpy import linalg as la
import sys
import math
import numpy.lib.scimath as nls
import scipy.constants as sci
sys.path.insert(0, '../../0/oliviabern')
from molecule import Molecule
sys.path.insert(0, '../../extra-files')
from masses import mass
from masses import charge

geom_string = open('../../extra-files/molecule.xyz').read()


h = open('../../extra-files/hessian.dat').read()
h = h.split()
for i in range(len(h)):
    h[i] = float(h[i])
h = np.array(h)
water = Molecule(geom_string)


def frequencies(water, h):
    w = str(water)
    N = int(w[0])
    h = h.reshape(3*N,3*N)
    print h
    atoms = []
    xyz = []
    lines = w.splitlines()
    for line in lines[2:]:
        atom, x, y, z = line.split()
        atoms.append(atom)
        xyz.append([float(x),float(y), float(z)])
    m = []
    molec = []
    for a in atoms:
        m.append([mass[charge[a]]]*3)
        molec.append(a)
    m = np.sqrt(m)
    M = np.diagflat(m)
    M = la.inv(M)
    H = np.dot(M,h) 
    H = np.dot(H,M)
    [K,Q] = la.eigh(H)
    q = np.dot(H,Q)
    q = np.transpose(q)
    k = K*(sci.physical_constants['hartree-joule relationship'][0])/((sci.physical_constants['atomic unit of length'][0])**2)/(sci.physical_constants['atomic mass unit-kilogram relationship'][0]) #converting to J, m2, and kg
    v = nls.sqrt(k)/(2*math.pi)
    V = v/(sci.physical_constants['speed of light in vacuum'][0]*100)
    a = 0
    b = 3
    out =''
    for l in range(len(V)):
        out += '{:d}\n'.format(N)
        out += '{:f}\n'.format(V[l])
        a = 0
        b = 3
        for i in range(N):
            out += '{:2s} {: >15.10f} {: >15.10f} {: >15.10f}'.format(molec[i], *xyz[i])
            out += '{: 15.10f} {: >15.10f} {: >15.10f}\n'.format(*q[l][range(a,b)])
            a += 3
            b += 3
        out += '\n' 
        
    file = open('output.xyz','w')
    file.write(out)
    file.close()
frequencies(water, h)
