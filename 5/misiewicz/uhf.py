# This project will be done in PsiAPI, rather than Psithon.

import collections
import itertools
import numpy as np
import psi4
import scipy.linalg

np.set_printoptions(threshold=np.inf, precision=5, linewidth=200, suppress=True)

def uhf(mol, tolerance=1.e-16):
    '''Run restricted Hartree-Fock on a molecule. Takes the molecule and the
    convergence tolerance as arguments. Returns the density matrix, the
    molecular energy, the coefficient matrix, the eigenvalues, the
    molecular integral tensor, and the number of unoccupied orbitals.'''
    
    if not isinstance(mol, psi4.core.Molecule):
        raise Exception("Hartree fock molecule must be a psi4 molecule.")
    if not isinstance(tolerance, float):
        raise Exception("Hartree-Fock tolerance must be a float.")
    
    psi4.set_options({"basis": "sto-3g"})
    
    # Construct the basis. Implicitly uses the basis option set earlier.
    basisset = psi4.core.BasisSet.build(mol)
    mints = psi4.core.MintsHelper(basisset)
    V_n = mol.nuclear_repulsion_energy()
    T = mints.ao_kinetic().to_array()
    T = scipy.linalg.block_diag(T, T)
    V = mints.ao_potential().to_array()
    V = scipy.linalg.block_diag(V, V)
    h = T + V
    S = mints.ao_overlap()
    S = S.to_array(S)
    S = scipy.linalg.block_diag(S, S)
    S = psi4.core.Matrix.from_array(S)
    S.power(-0.5, tolerance) # ...and now it's X. Remember, power returns a dimension!
    X = S.to_array()
    density_matrix = np.zeros(X.shape)
    I = mints.ao_eri().to_array()
    kron_matrix = np.zeros((2, 2, 2, 2))
    # TODO: There has to be a more elegant way to do this...
    kron_matrix[0, 0, 0, 0] =  kron_matrix[1, 1, 0, 0] =  kron_matrix[0, 0, 1, 1] =  kron_matrix[1, 1, 1, 1] = 1
    G = np.kron(kron_matrix, I)
    G = G.transpose((0, 2, 1, 3))
    
    # This code must run AFTER the basis set is built. Otherwise, the
    # atoms are marked as dummies and not included in natom().
    nocc = sum([int(mol.Z(n)) for n in range(mol.natom())]) - mol.molecular_charge()

    def iteration(D):
        v = np.einsum('abcd, db -> ac', G, D) - (
            np.einsum('abdc, db -> ac', G, D))
        f = h + v
        f_hat = X @ f @ X
        # Note that indexing is lost after this step - the eigenvectors are
        # put from lowest to highest eigenvalue, in C_hat.
        evals, C_hat = np.linalg.eigh(f_hat)
        # Define C as a matrix, so we have a conjugate transpose.
        C = np.asmatrix(X @ C_hat)
        C_full = C.copy()
        C[:, nocc:] = 0 # Zero terms from unoccupied orbitals.
        D = C @ C.getH() # See pg. 138 Szabo, 1st ed., for interpretation.
        # The columns of C no longer match our index, but the rows do.
        # And since D was defined in terms of the rows, D is properly indexed!
        E_e = np.trace((h + 1/2 * v) @ D)
        E = E_e + V_n
        return D, E, C_full, evals

    energy = 0
    energy_diff = tolerance + 1 # A dummy diff so the while suite executes.

    while energy_diff >= tolerance:
        previous_energy = energy
        density_matrix, energy, C, evals = iteration(density_matrix)
        energy_diff = abs(previous_energy - energy)

    return collections.namedtuple('hf_data', 'density,energy,C,evals,I,nocc')(
        density_matrix, energy, C, evals, G, nocc)

if __name__ == "__main__":
    psi4.set_memory('500 MB')
    h2o = psi4.geometry("""
1 2
O
H 1 0.96
H 1 0.96 2 104.5
""")
    hf_data = uhf(h2o)
    print(hf_data.energy)