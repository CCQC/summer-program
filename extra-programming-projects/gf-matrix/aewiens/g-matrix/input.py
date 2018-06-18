# This is the main program, which completes each part of the assignment using
# the modules "molecule" and "vib".
import numpy as np
from vib import IntCoVibAnalysis
from molecule import Molecule

# build HOOH and DOOD Molecule objects
hooh = Molecule(
"""
H
O 1 0.9625
O 2 1.4535 1 99.64
H 3 0.9625 2 99.64 1 113.7
"""
)

dood = Molecule(
"""
D
O 1 0.9625
O 2 1.4535 1 99.64
D 3 0.9625 2 99.64 1 113.7
"""
)

# build peroxide F matrix in C2-symmetrized internal coordinate basis
F = np.array(
#       s1,    s2,    s3,    s4,    s5,    s6
  [[ 8.201,-0.128,-0.021, 0.023,     0,     0], # s1
   [-0.128, 4.660, 0.825,-0.011,     0,     0], # s2
   [-0.021, 0.825, 1.035, 0.048,     0,     0], # s3
   [ 0.023,-0.011, 0.048, 0.031,     0,     0], # s4
   [     0,     0,     0,     0, 8.207,-0.036], # s5
   [     0,     0,     0,     0,-0.036, 0.869]] # s6
)

# build transformation matrix from regular internal coordinates (in zmat order)
# to C2-symmetrized internal coordinates for HOOH/DOD
n = 1 / np.sqrt(2)
U = np.array(
#      rOH,   rOO,  fHOO,  rOH', fHOO', tHOOH
  [[    +n,     0,     0,    +n,     0,     0], # s1
   [     0,    +1,     0,     0,     0,     0], # s2
   [     0,     0,    +n,     0,    +n,     0], # s3
   [     0,     0,     0,     0,     0,    +1], # s4
   [    +n,     0,     0,    -n,     0,     0], # s5
   [     0,     0,    +n,     0,    -n,     0]] # s6
)

# compute the G matrices for HOOH and DOOD

G = hooh.compute_g_matrix() # compute G matrix in unsymmetrized internal coordinate basis
m = hooh.masses
r = hooh.distance
f = hooh.angle
t = hooh.torsion
G = np.dot(U, np.dot(G, U.T)) # G -> U * G * U.T, transform to symmetrized internal coordinate basis
vib_hooh = IntCoVibAnalysis(F, G, blockdims=(4, 2))
print("HOOH G matrix")
print G.round(14)

G = dood.compute_g_matrix() # compute G matrix in unsymmetrized internal coordinate basis
G = np.dot(U, np.dot(G, U.T)) # G -> U * G * U.T, transform to symmetrized internal coordinate basis
vib_dood = IntCoVibAnalysis(F, G, blockdims=(4, 2))
print("DOOD G matrix")
print G.round(14)

# perform normal mode analysis and compute the harmonic frequencies

vib_hooh.print_cm_frequencies()
vib_dood.print_cm_frequencies()

# compute the total energy distributions of the normal modes

vib_hooh.print_total_energy_distribution()
vib_dood.print_total_energy_distribution()
