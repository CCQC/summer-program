#Project 1!
from p0 import Molecule
import numpy as np
import math
import masses
np.set_printoptions(suppress=True)

with open("hessian.dat", "r") as f:
    hess = f.readlines()

with open("molecule.xyz", "r") as f:
    coordinates = [line.rstrip() for line in f]


#1 Reading in the hessian
H = []
for i in hess:
    H.append(list(map(float, i.split())))
H = np.matrix(H)


#2 Building my molecule with the previously-made Molecule class

mol = Molecule("molecule.xyz")

def freq(mol, H):    
    
    #3 Mass-Weighted Hessian math (this was hard) 
    mass_lst = mol.mass_list
    N = len(mass_lst)
    M = np.zeros( (3*N,3*N), dtype = float )
    mw = [ mass_lst[int(i)] for i in range(N) for j in range(3) ]
    for i in range(len(mw)):
        M[i,i] = mw[i] ** -0.5
    
    F = M * H * M
    
    #4 Eigenvector and Eigenvalue time. Thanks numpy. 
    
    Val, Vec = np.linalg.eigh(F)
    
    
    #5 I do believe this will unmass the evecs!! 
    UwH = np.dot(M, Vec)
    
    #6 Figure out how to multiply these equations 
    
    #Let's make some easy unit conversion variables: 
    
    HtJ = 4.3597482e-18
    AMUtkg = 1.6603145e-27
    c = 29979245800 #cm/s (cause we want frequency!)
    BRtm = 5.29177e-11
    
    
    #k time
    
    k_1 = np.sqrt((abs(Val)*HtJ)/((BRtm**2)*AMUtkg))/(2*math.pi) 
    k_2 = k_1/c
    
    
    #Output!
    
    output = ""
    output_format = '{:2s}' + '{:  >15.10f}'*6 + '\n'
    for i in range(3*N):
        output +='{:d}\n{: >7.2f}  cm^-1\n'.format(N, np.absolute(k_2[i]))
        for j in range(N):
            atom = mol.atom_labels[j]
            x, y, z = mol.geom[j]
            x_1, y_1, z_1 = UwH[3*j:3*j+3, i]
            x_1, y_1, z_1 = float(x_1), float(y_1), float(z_1)
            output += '{:<2s}{:>20.12f}{:>20.12f}{:>20.12f}{:>20.12f}{:>20.12f}{:>20.12f}'.format(atom, x, y, z, x_1, y_1, z_1)
            output += "\n"
        output += "\n"
    
    filetime = open("FreqOutput.xyz", 'w')
    filetime.write(output)
    filetime.close
    
freq(mol, H)    
    
    
    
