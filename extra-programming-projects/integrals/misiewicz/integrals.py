import functools
import itertools
import math
import numpy as np
import os
import re
from scipy import special, misc
import sys
# TODO: Change summer-program into a package, so I don't have to do this.
sys.path.insert(0, os.path.join(sys.path[0], "..", "..", "..", "0", "misiewicz"))
import molecule
sys.path.insert(0, os.path.join(sys.path[0], "..", "..", "..", "extra-files")) 
import masses

# Uncomment when debugging.
np.set_printoptions(threshold=np.inf, precision=5, linewidth=200, suppress=True)

# Tuples nested three deep. First tuple separates atoms.
# Second tuple separates shells. Third tuple separates Gaussians.
# TODO: Replace this with a generator from a basis set file.
exponential_tuple = (
    (
        (3.42525091, 0.62391373, 0.16885540),
    ),
    (
        (6.36242139, 1.15892300, 0.31364979),
    ),
    (
        (16.1195750, 2.9362007, 0.7946505),
        (0.6362897, 0.1478601, 0.0480887)
    ),
    (
        (30.167871, 5.4951153, 1.4871927),
        (1.3148331, 0.3055389, 0.0993707)
    ),
    (
        (48.791113, 8.8873622, 2.405267),
        (2.2369561, 0.5198205, 0.1690618)
    ),
    (
        (71.616837, 13.045096, 3.5305122),
        (2.9412494, 0.6834831, 0.2222899)
    ),
    (
        (99.106169, 18.052312, 4.8856602),
        (3.7804559, 0.8784966, 0.2857144)
    ),
    (
        (130.70932, 23.808861, 6.4436083),
        (5.0331513, 1.1695961, 0.380389)
    ),
    (
        (166.67913, 30.360812, 8.2168207),
        (6.4648032, 1.5022812, 0.4885885)
    ),
    (
        (207.01561, 37.708151, 10.205297),
        (8.2463151, 1.9162662, 0.6232293)
    )
)

# The first index selects your subshell, while the second index selects your Gaussian.
# TODO: Replace this with a generator from a basis set file.
coefficient_tuple = (
    (0.15432897, 0.53532814, 0.44463454),
    (-0.09996723, 0.39951283, 0.70011547),
    (0.15591627, 0.60768372, 0.39195739)
)

# TODO: Replace this with a generator.
subshells = {
    "H": ["1s"],
    "He": ["1s"],
    "Li": ["1s", "2s"],
    "Be": ["1s", "2s"],
    "B": ["1s", "2s", "2p"],
    "C": ["1s", "2s", "2p"],
    "N": ["1s", "2s", "2p"],
    "O": ["1s", "2s", "2p"],
    "F": ["1s", "2s", "2p"],
    "Ne": ["1s", "2s", "2p"]
}

def subshell_to_orbitals(subshell):
    '''Takes a subshell as input and returns a tuple of all orbitals in the
    subshell as output.'''
    if subshell.endswith("s"):
        return (subshell,)
    if subshell.endswith("p"):
        shell = subshell[0]
        orbitals = ("pz", "px", "py") # Ordered to match psi4.
        return (shell + i for i in orbitals)
    return Exception("Invalid subshell.")

"""
Maybe this code would be useful for the generator? The heavy redundancies in
the basis set format make this difficult to write.

def read_nwchem_basisset(basis_filename):
    #Reads a basis set in NWChem format and converts it to ???.
    subshell = "\d[spdf]"
    mult_subshell = "((?:{},)*{})".format(subshell, subshell)
    format_dict = {"MS": mult_subshell}
    with open(basis_filename) as f:
        lines = return f.read().splitlines()
    for line in lines:
        # Check if the line defines ???. Capturing groups are in the
        # parentheses and brackets.
        match_obj = re.match("#BASIS SET: \({}\) -> \[{}\]".format(mult_subshell, mult_subshell)
        if match_obj:
            inp, out = match_obj.groups()
"""

def read_xyz(xyz_filename):
    """Input an xyz filename. Assumes units are in Angstroms.
    Outputs the atomic labels and the geometry matrix, in Bohr radii."""
    mol = molecule.Molecule(xyz_filename)
    mol.to_bohr()
    return mol.labels, np.asmatrix(mol.geom)

def generate_sto3g_data(atoms, geom):
    """Create the sto3g_data."""
    sto3g_data = []
    for atom, geom in zip(atoms, geom):
        for subshell, coefficient in zip(subshells[atom], coefficient_tuple):
            atom_index = masses.get_charge(atom) - 1
            for orbital in subshell_to_orbitals(subshell):
                # A full function will be needed for this if I extend into d
                # block elements.
                subshell_index = int(orbital[0]) - 1
                sto3g_data.append({
                    "a": exponential_tuple[atom_index][subshell_index],
                    "d": coefficient,
                    "r": geom.tolist()[0],
                    "l": 1 if orbital == "2px" else 0,
                    "m": 1 if orbital == "2py" else 0,
                    "n": 1 if orbital == "2pz" else 0,
                })
    return sto3g_data

def fact_double(x):
    """Inputs x and outputs (2x-1)!!."""
    return misc.factorial2(2*x-1, exact=True)

def fact(x):
    """Inputs x and outputs x!."""
    return misc.factorial(x, exact=True)

def f_j(j, l, m, a, b):
    """For two atoms and a directon, input an index over the sum of their
    angular momentum quantum numbers in that direction, their angular momentum
    quantum number in that direction, and the central displacement vector in
    that direction. Output a number. Used in matrix element calculations."""
    return sum([special.binom(l,k) * special.binom(m,j-k) * a ** (l-k) * b ** (m+k-j) for k in range(max(0, j-m), min(j,l)+1)])

def boys(nu, x):
    """Input nu and x and return the Boys function of nu and x."""
    if x >= 1.e-6:
        return 0.5 * x ** (-nu -0.5) * special.gamma(nu + 0.5) * special.gammainc(nu + 0.5, x)
    else:
        return 1 / (2 * nu + 1) - x / (2 * nu + 3)

def normalize(gaussian_list):
    """Given a list of gaussians, return the product of their normalization and 
    contraction coefficients."""
    def normalizer(a, l, m, n):
        j = l + m + n
        fact_term = np.prod([fact_double(j) for j in {l, m, n}])
        return (2 ** (2 * j) * a ** j / fact_term * (2 * a / np.pi) ** 1.5 ) ** 0.5

    norm = np.prod([normalizer(
        g["a"], g["l"], g["m"], g["n"]) * g["d"] for g in gaussian_list])

    return norm

def boys_manager(gaussian_list, func, boys_x, coord=[0, 0, 0]):
    """Given a list of Gaussians, a function to run on each Gaussian, a Boys
    function x value, and an optional coordinate, return the term that
    the Gaussians contribute to the matrix element."""
    # This needs to be a list to persist after iteration.
    gaussian_pairs = list(zip(gaussian_list[0::2], gaussian_list[1::2]))
    j_combo = {}

    # First, generate all possible (l, r) pairs.
    for p, n in (("l", 0), ("m", 1), ("n", 2)):

        lr_combos = []
        j_combo[p] = {}
        func_args = []

        for (g, h) in gaussian_pairs:
            # List all the Gaussian-dependent arguments.
            gamma = g["a"]+h["a"]
            p_v = (g["a"] * g["r"] + h["a"] * h["r"]) / (g["a"] + h["a"])
            func_args.extend([g[p], h[p], p_v[n] - g["r"][n], p_v[n] - h["r"][n], p_v[n] - coord[n], gamma])
            
            # Generate all possible (l, r) pairs for the ith pair.
            lr_combos.append(set())
            for l in range(g[p]+h[p]+1):
                for r in range(int(l/2)+1):
                    lr_combos[-1].add((l, r))

        # For each (l, r, l', r', ...) pair from the (l, r) pairs...
        # In practice, we have at most four tuple elements.
        for x in itertools.product(*lr_combos):
            lr_flattened = list(sum(x, ()))
            summa = sum(lr_flattened[::2])-2*sum(lr_flattened[1::2])
            for i in range(int(summa/2)+1):
                # ...generate all the +i pairs.
                index = tuple(lr_flattened + [i])
                # And hash the tuple to the function value.
                j_combo[p][index] = func(*index, *func_args)

    t = 0
    for l, m, n in itertools.product(j_combo["l"], j_combo["m"], j_combo["n"]):
        prod = j_combo["l"][l] * j_combo["m"][m] * j_combo["n"][n]
        l, m, n = list(l), list(m), list(n)
        i, j, k = l.pop(), m.pop(), n.pop()
        nu = sum(l[0::2]) + sum(m[0::2]) + sum(n[0::2]) - 2 * (sum(l[1::2]) + sum(m[1::2]) + sum(n[1::2])) - (i + j + k)
        t += boys(nu, boys_x) * prod

    return t

def gaussian_overlap(g, shift=(0, 0, 0)):
    """Compute and return the overlap integrals for Gaussians.
    The shift parameter shifts the effective angular momentum components."""

    # Negative angular momentum components mean the term with this integral
    # doesn't exist. Return 0.
    if any((i < 0 for i in {g[0]["l"], g[0]["m"], g[0]["n"], g[1]["l"], g[1]["l"], g[1]["l"]})):
        return 0

    gamma = g[0]["a"] + g[1]["a"]

    P = (g[0]["a"] * g[0]["r"] + g[1]["a"] * g[1]["r"]) / (g[0]["a"] + g[1]["a"])
    PA = P - g[0]["r"]
    PB = P - g[1]["r"]
    AB_2 = np.sum((g[0]["r"] - g[1]["r"]) ** 2)
    
    def coefficient_j(j):
        '''Assumes gamma is in scope.'''
        return (fact_double(j) / (2*gamma) ** j)
    
    def summation(ang_1, ang_2, index):
        '''Assumes orbital_1 and orbital_2 are in scope, as are PA and PB.'''
        PA_i, PB_i = PA[index], PB[index]
        return sum([f_j(2*j, ang_1, ang_2, PA_i, PB_i) * coefficient_j(j) for j in range(0, int((ang_1+ang_2)/2)+1)])
    
    def S_Null(ang_1, ang_2, index):
        '''Assumes gamma, PA, and PB are in scope.'''
        ang_2 += shift[index]
        return np.sqrt(np.pi/gamma) * summation(ang_1, ang_2, index)
    
    prod = np.prod([S_Null(g[0][j], g[1][j], i) for j, i in {("l", 0), ("m", 1), ("n", 2)}])
    
    return np.exp(- g[0]["a"] * g[1]["a"] * AB_2 / gamma) * prod

def gaussian_recursion(orbital_list, gaussian_list, func, S, index):
    """From a list of orbitals, add the specified function of each combination
    of Gaussians (one per orbital) to S[index]. gaussian_list should receive an
    empty list when called by a function other than itself."""
    # WARNING: If you possibly can, never use a recursive for loop.
    if orbital_list: # If we have orbitals to represent in the function...
        o = orbital_list.pop(0) # Pick the first orbital...
        for (a, d) in zip(o["a"], o["d"]): # For each Gaussian...
            gaussian_list.append({
                "a": a, "l": o["l"], "m": o["m"], "n": o["n"], "d": d,
                "r": np.asarray(o["r"])})
            # list() makes the orbital_list arg a new object in the next layer,
            # so it doesn't interfere with this list. Otherwise, the next time
            # this loops, orbital_list is depleted.
            # gaussian_list is filled by recursion, so it needs no new object.
            gaussian_recursion(list(orbital_list), gaussian_list, func, S, index)
            gaussian_list.pop()
    else:
        S[index] += (normalize(gaussian_list) * func(gaussian_list))

def kinetic_term(gaussian_list):
    """Given a list of Gaussians for a particular inde of the electronic
    kinetic matrix, return their contribution to the value at the index."""
    # Slide 14 at http://folk.uio.no/helgaker/talks/SostrupIntegrals_10.pdf
    # suggests a simpler way to do this calculation. Confirm?
    
    g1, g2 = gaussian_list
    
    def k_gaussian(shift):
        return gaussian_overlap(gaussian_list, shift)
    
    t = g2["a"] * (2 * (g2["l"] + g2["m"] + g2["n"]) + 3) * k_gaussian((0, 0, 0))
    t -= 2 * g2["a"] ** 2 * k_gaussian((2, 0, 0))
    t -= 2 * g2["a"] ** 2 * k_gaussian((0, 2, 0))
    t -= 2 * g2["a"] ** 2 * k_gaussian((0, 0, 2))
    t -= 0.5 * g2["l"] * (g2["l"] - 1) * k_gaussian((-2, 0, 0))
    t -= 0.5 * g2["m"] * (g2["m"] - 1) * k_gaussian((0, -2, 0))
    t -= 0.5 * g2["n"] * (g2["n"] - 1) * k_gaussian((0, 0, -2))
    
    return t

def attraction_term(Z, C, g):
    """Given a list of Gaussians for a particular index of the nuclear attraction matrix, return their contribution to the value at the index."""
    gamma = g[0]["a"] + g[1]["a"]
    AB_2 = np.sum((g[0]["r"] - g[1]["r"]) ** 2)
    P = (g[0]["a"] * g[0]["r"] + g[1]["a"] * g[1]["r"]) / (g[0]["a"] + g[1]["a"])
    PC = P - C
    PC_2 = np.sum(PC ** 2)
    
    coefficient = Z * -2 * np.pi / gamma * np.exp(- g[0]["a"] * g[1]["a"] / gamma * AB_2)
    
    def v(l, r, i, la, lb, pa, pb, pc, gamma):
        return (-1) ** l * f_j(l, la, lb, pa, pb) * (-1) ** i * fact(l) * pc ** (l - 2 * r - 2 * i) * (1 / 4 / gamma) ** (r + i) / fact(r) / fact(i) / fact(l - 2 * r - 2 * i)

    return coefficient * boys_manager(g, v, gamma*PC_2, C)

def repulsion_term(g):
    """Given a list of Gaussians for a particular index of the electron repulsion matrix, return their contribution to the value at the index."""
    gam_1, gam_2 = g[0]["a"] + g[1]["a"], g[2]["a"] + g[3]["a"]
    AB_2, CD_2 = np.sum((g[0]["r"] - g[1]["r"]) ** 2), np.sum((g[2]["r"] - g[3]["r"]) ** 2)
    P = (g[0]["a"] * g[0]["r"] + g[1]["a"] * g[1]["r"]) / (g[0]["a"] + g[1]["a"])
    Q = (g[2]["a"] * g[2]["r"] + g[3]["a"] * g[3]["r"]) / (g[2]["a"] + g[3]["a"])
    #PA, PB, QC, QD = P - g[0]["r"], P - g[1]["r"], Q - g[2]["r"], Q - g[3]["r"]
    PQ_2 = np.sum((P - Q) ** 2)
    delta = 1/4 * (1/gam_1 + 1/gam_2)
    
    coefficient = 2 * np.pi ** 2 / gam_1 / gam_2 * (np.pi / (gam_1 + gam_2)) ** 0.5 * np.exp(- g[0]["a"] * g[1]["a"] / gam_1 * AB_2) * np.exp(- g[2]["a"] * g[3]["a"] / gam_2 * CD_2)
    
    def g_prod(l, r, lp, rp, i, la, lb, pa, pb, px, gam_1, lc, ld, qc, qd, qx, gam_2):
        def theta(l, la, lb, a, b, r, gamma):
            return f_j(l, la, lb, a, b) * fact(l) * gamma ** (r-l) / fact(r) / fact(l-2*r)
        return (-1) ** l * theta(l, la, lb, pa, pb, r, gam_1) * theta(lp, lc, ld, qc, qd, rp, gam_2) * (-1) ** i * (2 * delta) ** (2 * (r + rp)) * fact(l + lp - 2 * r - 2 * rp) * delta ** i * (px-qx) ** (l + lp - 2 * (r + rp + i)) / (4 * delta) ** (l + lp) / fact(i) / fact(l + lp - 2 * (r + rp +i))
    
    boys_x = PQ_2 / 4 / delta

    return coefficient * boys_manager(g, g_prod, boys_x)

def overlap_matrix(sto3g_data):
    """Compute the orbital overlap matrix."""
    number_ao = len(sto3g_data)
    S = np.zeros((number_ao, number_ao))
    for index in np.ndindex(S.shape):
        orbitals = [sto3g_data[i] for i in index]
        gaussian_recursion(orbitals, [], gaussian_overlap, S, index)

    return S

def kinetic_matrix(sto3g_data):
    """Compute the electron kinetic energy matrix."""
    number_ao = len(sto3g_data)
    T = np.zeros((number_ao, number_ao))
    for index in np.ndindex(T.shape):
        orbitals = [sto3g_data[i] for i in index]
        gaussian_recursion(orbitals, [], kinetic_term, T, index)

    return T

def nuclear_attraction_matrix(sto3g_data, labels, geom):
    """Compute the nuclear attraction matrix."""
    number_ao = len(sto3g_data)
    V = np.zeros((number_ao, number_ao))
    for label, geom in zip(labels, geom):
        Z = masses.get_charge(label)
        V_temp = np.zeros((number_ao, number_ao))
        for index in np.ndindex(V.shape):
            orbitals = [sto3g_data[i] for i in index]
            n_attract = functools.partial(attraction_term, Z, geom.tolist()[0])
            gaussian_recursion(orbitals, [], n_attract, V_temp, index)

        V += V_temp
        
    return V

def repulsion_tensor(sto3g_data):
    """Compute the repulsion tensor."""
    number_ao = len(sto3g_data)
    G = np.zeros((number_ao, number_ao, number_ao, number_ao))
    for index in np.ndindex(G.shape):
        orbitals = [sto3g_data[i] for i in index]
        gaussian_recursion(orbitals, [], repulsion_term, G, index)
    return G

def save_repulsion_tensor(filename, repulsion_tensor):
    """Save the repulsion tensor to the specified file."""
    with open("S."+filename+".npy", "bw") as f:
        np.save(f, repulsion_tensor)

if __name__ == "__main__":
    
    labels, geom = read_xyz("h2o.xyz")
    sto3g_data = generate_sto3g_data(labels, geom)
    overlap = overlap_matrix(sto3g_data)
    kinetic = kinetic_matrix(sto3g_data)
    nuclear = nuclear_attraction_matrix(sto3g_data, labels, geom)
    repulsion = repulsion_tensor(sto3g_data)
    save_repulsion_tensor("h2o", repulsion)
