import sys
sys.path.insert(0, '../../extra-files/')
sys.path.insert(0, '../../0/huples-cat/')
from masses import get_charge, get_mass
from molecule import Molecule
import numpy as np

#open hessian
lines = open('../../extra-files/hessian.dat').readlines()

#create molecule object
mol = Molecule('../../extra-files/molecule.txt')

#create hessian numpy array
H = []
for line in lines:
    H.append = line.split()
H = np.array(H, dtype = float)

def frequency(Molecule = mol, np.array = H):
    
