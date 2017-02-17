from psi4.core import BasisSet, MintsHelper
import numpy as np
from scipy import linalg as la
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
    return np.einsum('Pp, Pqrs->pqrs', C,
                np.einsum('Qq,PQrs->Pqrs', C,
                    np.einsum('Rr,PQRs->PQrs', C,
                        np.einsum('Ss,PQRS->PQRs', C, gao))))

def transform_integrals_einsum_noddy(C, gao):
    """
    Einsum version of integral transformation
    Noddy algorithm
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


def transform_integrals_df(C, scf, df_basis):
    """
    Generate density-fitted integrals
    """
    # TODO: Remove this line
    np.set_printoptions(precision=3, linewidth=150)

    # Form (pq|P)
    pqP = eri(scf.basis, scf.basis, df_basis)
    pqP = block_tei_3c(pqP)
    #print('pqP')
    #print(pqP)

    # Form X≡ J^{-1/2}
    J = eri(df_basis, df_basis)
    X = la.inv(la.sqrtm(J))
    #X = la.block_diag(X, X)
    #print('X')
    #print(X)

    #print('C')
    #print(C)
    #print(C.shape, pqP.shape)
    b_pqP = np.einsum('pqQ,PQ->pqP', pqP, X)
    b_pnP = np.einsum('pm,mnP->pnP', C, b_pqP)
    b_pqP = np.einsum('qn,pnP->pqP', C, b_pnP)
    gmo = np.einsum('pqP, rsP->pqrs', b_pqP, b_pqP)

    return gmo.transpose(0, 2, 1, 3) - gmo.transpose(0, 2, 3, 1)


def block_tei(a):
    """
    Spin block two-electron integrals
    """
    I = np.eye(2)
    a = np.kron(I, a)
    return np.kron(I, a.T)


def block_tei_3c(a):
    """
    Spin block 3-center two-electron integrals
    """
    I = np.eye(2)
    return np.kron(I, a.T).T


def eri(bs1, bs2, bs3=0, bs4=0):
    zero = BasisSet.zero_ao_basis_set()
    mints = MintsHelper(bs1)

    if not bs3:
        res = mints.ao_eri(bs1, zero, bs2, zero)
    elif not bs4:
        res = mints.ao_eri(bs1, bs2, bs3, zero)

    return np.squeeze(res)
