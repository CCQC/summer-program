import numpy as np
import psi4
import sys
sys.path.append('../../3/mdav2')
import rhf

def mp2(C, g, e, nocc):
#Calculates mp2 (moller-plesset @ 2nd order) correlation energy.
#Takes in a coefficient matrix (C), two-electron integrals (g),
#orbital energies (e), and an occupation number (nocc)
    phi = np.einsum('pQRS,pP->PQRS',
          np.einsum('pqRS,qQ->pQRS',
          np.einsum('pqrS,rR->pqRS',
          np.einsum('pqrs,sS->pqrS', g, C), C), C), C)
#translate to physicists notation
    phi_alt = phi.transpose(0,1,3,2)
    mp2_correction = 0
#Calculate mp2 corr energy via the RMP2 simplifications
    for I in range(nocc):
        for J in range(nocc):
            for A in range(nocc,len(e)):
                for B in range(nocc, len(e)):
                    mp2_correction += (phi[I][J][A][B]*(2*phi[I][J][A][B]
                                   - phi_alt[I][J][A][B]))\
                                   /(e[I]+e[J]- e[A] - e[B])
    return mp2_correction
    
if __name__ == '__main__':
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    EHF, msg, C, g, e, nocc = rhf.rhf(mol, basis='cc-sto-3g', iterlim=500)
    mp2_corr = mp2(C, g, e, nocc)
    E =EHF + mp2_corr
    print("EHF: {}\nMP2 CORRECTION: {}\nE=(EHF+MP2): {}".format(EHF, mp2_corr, E))
 
