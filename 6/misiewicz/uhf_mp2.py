import itertools
import numpy as np
import psi4
import sys
sys.path.insert(0, '../../5/misiewicz')
import uhf
import utility

def mp2(C, evals, int_tensor, nocc):
    dim = C.shape[0]
    nvir = dim - nocc
    # We only need tensor elts of form <occ, occ, vir, vir>. Construct them.
    C_occ = C[:,:nocc].conjugate() # Grab the columns for occupied orbitals.
    C_vir = C[:,nocc:] # Grab the columns for virtual orbitals.
    aint = utility.tensor_basis_transform(int_tensor, C_occ, C_vir)
    e = 0
    for i, j, a, b in itertools.product(range(nocc), range(nocc), range(nvir), range(nvir)):
        double_int = aint[i,j,a,b] - aint[i,j,b,a]
        e += double_int * np.conjugate(double_int) / (
        	evals[i] + evals[j] - evals[nocc+a] - evals[nocc+b])
    
    return e / 4

if __name__ == "__main__":
    psi4.set_memory('500 MB')
    h2o = psi4.geometry("""
1 2
O
H 1 0.96
H 1 0.96 2 104.5
""")
    hf_data = uhf.uhf(h2o)[2:] # Throw away first two data points...
    print(mp2(*hf_data))