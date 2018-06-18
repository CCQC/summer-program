from molecule import Molecule

hooh = Molecule(
"""
H
O 1 0.9625
O 2 1.4535 1 99.64
H 3 0.9625 2 99.64 1 113.7
"""
)

print( hooh.compute_g_matrix() )

