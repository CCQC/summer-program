import numpy
import geom
import math

def readhessian(filename):
    hessianfile = open(filename, 'r')
    hessianlines = hessianfile.readlines()
    
    
    hessianmatrix = numpy.zeros(shape=(len(hessianlines),len(hessianlines)))
    for line in range(0,len(hessianlines)):
        hessianmatrix[line] = hessianlines[line].split() 
    return hessianmatrix 
def frequencies(molecule,hessian):
    
    weightedhessian = numpy.zeros(shape=hessian.shape)
    
    for A in range(0,molecule.natom * 3):
        for B in range(0, molecule.natom * 3):
            weightedhessian[A][B] = hessian[A][B]/math.sqrt(molecule.masses[math.floor(A/3.)]*molecule.masses[math.floor(B/3.)])
    print(weightedhessian) 
    

    

        
    
