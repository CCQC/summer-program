import numpy as np
from itertools import product


def transform_integrals_noddy(C, gao):
    """
    Transform the integrals from the AO basis to the MO basis
    """
    nbf = len(gao)
    gmo = np.zeros(gao.shape)
    # Noddy algorithm, O(N^8)
    for p in range(nbf):
        for q in range(nbf):
            for r in range(nbf):
                for s in range(nbf):
                    for mu in range(nbf):
                        for nu in range(nbf):
                            for rho in range(nbf):
                                for sigma in range(nbf):
                                    gmo[p, q, r, s] += gao[mu, nu, rho, sigma]*\
                                            C[mu, p]*C[nu, q]*C[rho, r]*C[sigma, s]
    return gmo

def transform_integrals_einsum(C, gao):
    """
    Einsum version of integral transformation
    """
    return np.einsum('PQRS,Pp,Qq,Rr,Ss->pqrs', gao, C, C, C, C)

def transform_integrals(C, gao):
    """
    Transform integrals the efficient way O(N^5)
    """
    nbf = len(gao)

    g = np.zeros(gao.shape)
    for mu in range(nbf):
        for nu in range(nbf):
            for rho in range(nbf):
                for s in range(nbf):
                    for sigma in range(nbf):
                        g[mu, nu, rho, s] += gao[mu, nu, rho, sigma]*C[sigma, s]

    gmo = np.zeros(gao.shape)
    for mu in range(nbf):
        for nu in range(nbf):
            for r in range(nbf):
                for s in range(nbf):
                    for rho in range(nbf):
                        gmo[mu, nu, r, s] += g[mu, nu, rho, s]*C[rho, r]

    g.fill(0)
    for mu in range(nbf):
        for q in range(nbf):
            for r in range(nbf):
                for s in range(nbf):
                    for nu in range(nbf):
                        g[mu, q, r, s] += gmo[mu, nu, r, s]*C[nu, q]

    gmo.fill(0)
    for p in range(nbf):
        for q in range(nbf):
            for r in range(nbf):
                for s in range(nbf):
                    for mu in range(nbf):
                        gmo[p, q, r, s] += g[mu, q, r, s]*C[mu, p]

    return gmo

def transform_integrals_itertools(C, gao):
    """
    Simplify iteration using a cartesian product
    """
    nbf = len(gao)
    g1, g2, g3, g4 = np.zeros(gao.shape), np.zeros(gao.shape), np.zeros(gao.shape), np.zeros(gao.shape)

    for μ, ν, ρ, σ, s in product(range(nbf), repeat=5):
        g1[μ, ν, ρ, s] += gao[μ, ν, ρ, σ]*C[σ, s]
    for μ, ν, ρ, r, s in product(range(nbf), repeat=5):
        g2[μ, ν, r, s] += g1[μ, ν, ρ, s]*C[ρ, r]
    for μ, ν, q, r, s in product(range(nbf), repeat=5):
        g3[μ, q, r, s] += g2[μ, ν, r, s]*C[ν, q]
    for μ, p, q, r, s in product(range(nbf), repeat=5):
        g4[p, q, r, s] += g3[μ, q, r, s]*C[μ, p]

    return g4
