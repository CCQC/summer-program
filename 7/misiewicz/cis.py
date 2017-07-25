import numpy as np
import sys
sys.path.insert(0, '../../5/misiewicz')
import uhf
sys.path.insert(0, '../../6/misiewicz')
import utility
import scipy.linalg as la

def index_splitter(index, nocc):
    return index // nocc + nocc, index % nocc

def cis(C, evals, int_tensor, nocc):
    norb = C.shape[0]
    nvir = norb - nocc
    aint = utility.tensor_basis_transform(int_tensor, C.conjugate(), C)
    H = np.zeros((nvir * nocc, nvir * nocc))
    for x, y in np.ndindex(H.shape):
        # Do all excitations of a given ground orbital, then move to the next ground orbital.
        a, i = index_splitter(x, nocc)
        b, j = index_splitter(y, nocc)
        val = aint[a, j, i, b] - aint[a, j, b, i]
        if (i, a) == (j, b):
            val += (evals[a] - evals[i])
        H[x, y] = val
    return np.linalg.eigh(H)

def print_cis_data(evals, C, nocc):
    # TODO: Clean up my string formatting.
    strings = ["State", " " * 3, "Exc. Energy", " " * 3, "Sig. Excitations"]
    print(("{}" * len(strings)).format(*strings) + "\n")
    format_str = "{:>" + str(len(strings[0])) + "}"  "{:>" + str(len(strings[1])) + "}" + "{:" +str(len(strings[2])) + ".8f}" + "{:>" + str(len(strings[3])) + "}"
    # For each excited state, take its number, eigenvalue, and eigenvector...
    for k, (eival, eivec) in enumerate(zip(evals, C.transpose()), start=1):
        k_str = format_str.format(k, strings[1], eival, strings[3])
        # Search for significant excitations
        for i, n in enumerate(eivec):
            if n**2 > 0.1:
                vir, occ = index_splitter(i, nocc)
                k_str += "{:>3} -> {:>3} @ {:.3f}%".format(occ, vir, n**2*100)
        print(k_str)
    	

if __name__ == "__main__":
    import psi4
    psi4.set_memory('500 MB')
    h2o = psi4.geometry("""
  O
  H 1 1.1
  H 1 1.1 2 104.0
""")
    hf_data = uhf.uhf(h2o)[2:] # Throw away first two data points...
    print_cis_data(*cis(*hf_data), hf_data[-1])
