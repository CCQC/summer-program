import numpy as np
from scipy import linalg as la
import sys

sys.path.append('../../')

import conv as cnv

from fcallen import F,G

GF = np.dot(G,F)
print(GF)
print(la.norm(GF))
		
np.savetxt("gf-matrix.dat",GF,"%15.7f"," ","\n")

e,l = la.eigh( GF )

for i in e:
	print(i)

"""
conv =  np.sqrt(cnv.hartree2J/(cnv.amu2kg*cnv.bohr2m**2))/(cnv.c*2*np.pi)

frequencies = [conv*np.sqrt(i) for i in e]

for k in frequencies:
	print(k)
"""
