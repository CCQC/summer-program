
import numpy as np
from vib import IntCoVibAnalysis
from molecule import Molecule

h2cn = Molecule(
"""
N
C 1 1.253338969095741
H 2 1.095130571925020 1 121.213471579853874
H 2 1.095130571925020 1 121.213471579853874 3 180.000000000000000
"""
)


G = h2cn.compute_g_matrix()  ## unsymmetrized internal coordinate basis

print(G)
