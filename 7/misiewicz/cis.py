import numpy as np
import sys
sys.path.insert(0, '../../5/misiewicz')
import uhf
sys.path.insert(0, '../../6/misiewicz')
import utility
import scipy.linalg as la

def cis(C, evals, int_tensor, nocc):
    norb = C.shape[0]
    nvir = norb - nocc
    aint = utility.tensor_basis_transform(int_tensor, C.conjugate(), C)
    H = np.zeros((nvir * nocc, nvir * nocc))
    for x, y in np.ndindex(H.shape):
        # Do all excitations of a given ground orbital, then move to the next ground orbital.
        a, i = x // nocc + nocc, x % nocc
        b, j = y // nocc + nocc, y % nocc
        val = aint[a, j, i, b] - aint[a, j, b, i]
        if (i, a) == (j, b):
            val += (evals[a] - evals[i])
        H[x, y] = val
    return np.linalg.eigh(H)
    
if __name__ == "__main__":
    import psi4
    psi4.set_memory('500 MB')
    h2o = psi4.geometry("""
  O
  H 1 1.1
  H 1 1.1 2 104.0
""")
    hf_data = uhf.uhf(h2o)[2:] # Throw away first two data points...
    evals, C = cis(*hf_data)
    print(evals)
