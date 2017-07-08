# This project will be done in PsiAPI, rather than Psithon.

import itertools
import numpy as np
import psi4
psi4.set_memory('500 MB')

h2o = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

psi4.set_options({"basis": "sto-3g"})
tolerance = 1.e-16

# Define the basis of h2o. Implicitly uses the basis option set earlier.
basisset = psi4.core.BasisSet.build(h2o)
mints = psi4.core.MintsHelper(basisset)
V_n = h2o.nuclear_repulsion_energy()
T = mints.ao_kinetic().to_array()
V = mints.ao_potential().to_array()
h = T + V
X = mints.ao_overlap() # This is actually S...
X.power(-0.5, 1.e-16) # ...and now it's X. Remember, power returns a dimension!
X = X.to_array()
D = np.zeros(X.shape)

# This code must run AFTER the basis set is built. Otherwise, the
# atoms are marked as dummies and not included in natom().
electron = sum([int(h2o.Z(n)) for n in range(h2o.natom())]) - h2o.molecular_charge()
if electron % 2:
	raise Exception("Must have a closed-shell system for RHF.")
else:
	nocc = electron // 2

def iteration(D, previous=0):
	# Convert from chemist to physicist notation.
	I = mints.ao_eri().to_array().transpose((0, 2, 1, 3))
	v = np.einsum('abcd, db -> ac', I, D) - (
		0.5 * np.einsum('abdc, db -> ac', I, D))
	f = h + v
	f_hat = X @ f @ X
	evals, C_hat = np.linalg.eigh(f_hat)
	# Define C as a matrix, so we have a conjugate transpose.
	C = np.asmatrix(X @ C_hat)
	# Zero terms from unoccupied orbitals.
	C[:, nocc:] = 0
	D = 2 * C @ C.getH()
	E_e = np.trace((h + 1/2 * v) @ D)
	E = E_e + V_n
	return D, E, abs(previous - E)

D, E, _ = iteration(D)

diff = tolerance + 1 # A dummy diff so the code executes at least once.
while diff >= tolerance:
    D, E, diff = iteration(D, E)

print(D)
print(E)