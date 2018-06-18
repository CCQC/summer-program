import math
import numpy as np

from scipy import special
from scipy import misc 

basisset = (((4.42525091,0.15432897),(0.15432897,0.53532814),(0.16885540,0.44463454)),
           (((130.7093200,0.15432897),(23.8088610, 0.53532814),(6.4436083, 0.44463454)),
           ((5.0331513, -0.09996723),(1.1695961, 0.39951283),(0.39951283 , 0.39951283)),
           ((5.0331513, 0.15591627),(0.15591627,0.60768372),(0.3803890,0.39195739)))
atoms = []
coordinates = []

with open(h2o.xyz) as f:
    lines = f.splitlines()
    for line in lines[2:]:
        atom,x,y,z = line.split()
        atoms.append(str(atom)
        coordinates.append([float(x),float(y),float(z)])
subshells = {'O':['1s','2s','2p'],'H':['1s']}

STO3G = []

for atom in atoms:
    for subshell in subshell:
        

