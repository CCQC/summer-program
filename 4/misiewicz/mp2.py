# http://www.helsinki.fi/kemia/fysikaalinen/opetus/jlk/luennot/Lecture4.pdf
# Technically, the 2*IJAB-IJBA integral should be complex-conjugated.
# Currently, these integrals are always real, but if somebody chose an
# imaginary basis for some bizarre reason, an imaginary integral is possible.
# I include the conjugation here for the sake of mathematical rigor.

import itertools
import numpy as np
import psi4
import sys
sys.path.insert(0, '../../3/misiewicz')
import hartree_fock
sys.path.insert(0, '../../6/misiewicz')
import utility

def mp2(C, evals, int_tensor, nocc):
    dim = C.shape[0]
    nvir = dim - nocc
    C_occ = C[:,:nocc].conjugate() # Grab the columns for occupied orbitals.
    C_vir = C[:,nocc:] # Grab the columns for virtual orbitals.
    aint = utility.tensor_basis_transform(int_tensor, C_occ, C_vir)
    e = 0
    for i, j, a, b in itertools.product(range(nocc), range(nocc), range(nvir), range(nvir)):
        e += aint[i,j,a,b] * (
        	np.conjugate(2 * aint[i,j,a,b] - aint[i,j,b,a])) / (
        	evals[i] + evals[j] - evals[nocc+a] - evals[nocc+b])
    
    return e

if __name__ == "__main__":
    psi4.set_memory('500 MB')
    h2o = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")
    hf_data = hartree_fock.rhf(h2o)[2:] # Throw away first two data points...
    print(mp2(*hf_data))