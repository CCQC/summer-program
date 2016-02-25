#!/usr/local/bin/python

#######################################################
### 4 helper functions from psi4 for RHF class ########
######################################################


import psi4

def get_nelec(mol):
####  get # of occupied SPIN orbitals from psi4 molecule class
    char = mol.molecular_charge()
    nelec = -char
    for A in range(mol.natom()):
        nelec += mol.Z(A)
    return int(nelec)

def get_docc(mol):
    return int( get_nelec(mol)/2 )

def get_nbf(mints):
# get # of basis functions from psi4 mints library
    nbf = mints.basisset().nbf()
    if( nbf**4 *8.0 > psi4.get_memory()  ):
        raise Exception("Integrals exceed memory limit.")
    return nbf

def get_conv():
# get convergence
    return psi4.get_global_option('E_CONVERGENCE')

def get_maxiter():
# get max # of iterations
    return psi4.get_option('SCF','MAXITER')
