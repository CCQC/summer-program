import numpy
import geom
import math

def readhessian(filename):
    """opens and reads a hessian matrix from the given filename"""
    hessianfile = open(filename, 'r')
    hessianlines = hessianfile.readlines()
    
    
    hessianmatrix = numpy.zeros(shape=(len(hessianlines),len(hessianlines)))
    for line in range(0,len(hessianlines)):
        hessianmatrix[line] = hessianlines[line].split() 
    
    return hessianmatrix 
def frequencies(molecule,hessian):
    #calculates normal modes and frequencies from a Hessian matrix. Takes a molecule file and a \
    # 3*natom X 3*natom hessian numpy array as input.
    weightedhessian = numpy.zeros(shape=hessian.shape)
    
    for A in range(0,molecule.natom * 3):
        for B in range(0, molecule.natom * 3):
            #print(A, B)
            #print(molecule.masses[A//3])
            #print(molecule.masses[math.floor(A/3.)])
            weightedhessian[A][B] = hessian[A][B]/math.sqrt(molecule.masses[math.floor(A/3.)]*molecule.masses[math.floor(B/3.)])
    
    evals,evecs = numpy.linalg.eig(weightedhessian)

    unweighted = numpy.zeros(shape=(molecule.natom * 3,molecule.natom * 3))

    for vector in range(0,molecule.natom * 3):
        unweighted[vector] = evecs[vector]/math.sqrt(molecule.masses[math.floor(vector/3.)])
    reorg = unweighted.T
    return reorg, evals

        
    
