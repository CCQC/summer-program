
from molecule import Molecule

h2o = Molecule(
"""
O
H 1 1.09520
H 1 1.09520 2 109.0000
"""
)

print( h2o.compute_g_matrix() )

