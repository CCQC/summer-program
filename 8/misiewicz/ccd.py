import itertools
import numpy as np
import psi4
import sys
sys.path.insert(0, '../../5/misiewicz')
import uhf
sys.path.insert(0, '../../6/misiewicz')
import utility

def ccd(C, evals, int_tensor, nocc):
    tolerance=1.e-10
    norb = C.shape[0]
    nvir = norb - nocc
    aint = utility.tensor_basis_transform(int_tensor, C.conjugate(), C)
    T = np.zeros((nocc, nocc, nvir, nvir))
    anti_tensor = aint - aint.transpose(0,1,3,2)
    energy = 0
    energy_diff = tolerance + 1 # A dummy diff so the while suite executes.

    o_sl = slice(None, nocc)
    v_sl = slice(nocc, None)

    def s(a):
        return a + nocc

    while energy_diff >= tolerance:
        previous_energy = energy
        T_new = np.zeros((nocc, nocc, nvir, nvir))
        # If nocc = nvir, Andreas's method could cause problems...
        t_sum_0 = anti_tensor.transpose(2, 3, 0, 1)[o_sl, o_sl, v_sl, v_sl] + 0.5 * np.einsum('ijcd, abcd -> ijab', T, anti_tensor[v_sl, v_sl, v_sl, v_sl]) + 0.5 * np.einsum('klab, klij -> ijab', T, anti_tensor[o_sl, o_sl, o_sl, o_sl]) + (
        np.einsum('jkbc, akic -> ijab', T, anti_tensor[v_sl, o_sl, o_sl, v_sl]) - np.einsum('ikbc, akjc -> ijab', T, anti_tensor[v_sl, o_sl, o_sl, v_sl]) - np.einsum('jkac, bkic -> ijab', T, anti_tensor[v_sl, o_sl, o_sl, v_sl]) + np.einsum('ikac, bkjc -> ijab', T, anti_tensor[v_sl, o_sl, o_sl, v_sl])
        ) - 0.5 * (
        np.einsum('ijac, klbd, klcd -> ijab', T, T, anti_tensor[o_sl, o_sl, v_sl, v_sl]) - np.einsum('ijbc, klad, klcd -> ijab', T, T, anti_tensor[o_sl, o_sl, v_sl, v_sl])
        ) - 0.5 * (
        np.einsum('ikab, jlcd, klcd -> ijab', T, T, anti_tensor[o_sl, o_sl, v_sl, v_sl]) - np.einsum('jkab, ilcd, klcd -> ijab', T, T, anti_tensor[o_sl, o_sl, v_sl, v_sl])
        ) + 0.25 * np.einsum('ijcd, klab, klcd -> ijab', T, T, anti_tensor[o_sl, o_sl, v_sl, v_sl]) + (
        np.einsum('ikac, jlbd, klcd -> ijab', T, T, anti_tensor[o_sl, o_sl, v_sl, v_sl]) - np.einsum('jkac, ilbd, klcd -> ijab', T, T, anti_tensor[o_sl, o_sl, v_sl, v_sl]))
        for i, j, a, b in np.ndindex(T_new.shape):
            T[i, j, a, b] = t_sum_0[i, j, a, b] / (evals[i] + evals[j] - evals[s(a)] - evals[s(b)])
        energy = np.sum(anti_tensor[o_sl, o_sl, v_sl, v_sl] * T)/4
        print(energy)
        energy_diff = abs(previous_energy - energy)

    return energy

if __name__ == "__main__":
    psi4.set_memory('500 MB')
    h2o = psi4.geometry("""
1 2
O
H 1 1.1
H 1 1.1 2 104.0
""")
    hf_data = uhf.uhf(h2o)[2:] # Throw away first two data points...
    # My answers differ from Andreas's in the tenths place.
    # Likely cause: He did more UHF iterations than I did.
    print(ccd(*hf_data))