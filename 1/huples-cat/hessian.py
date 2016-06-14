import sys
sys.path.insert(0, '../../extra-files/')
sys.path.insert(0, '../../0/huples-cat/')
from masses import get_charge, get_mass
from molecule import Molecule

#open Hessian
lines = open('../../extra-files/hessian.dat').readlines()

#create molecule object
mol = Molecule('../../extra-files/molecule.txt')

