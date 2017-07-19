import numpy
import geom
import math

ccms = 3E10 #speed of light in centimeters per second
eh = 4.36E-18 #J/hartree
ao = 5.29E-11 #m/br
u = 1.66E-27 #kg/amu

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
    evals = evals * (eh/(u * ao **2))
    frequencies = []
    for idx,frequency in enumerate(evals):
        freq = ( (1./(math.pi * 2) * numpy.sqrt(frequency + 0J))/ccms) #+0J to force complex evaluation
        #turns imaginary frequencies to negative
        if freq.imag != 0:
            freq = -1*freq.imag

        frequencies.append(freq)
    sortingdictionary = dict(zip(frequencies,reorg))
    modes = []
    freqs = []
    #sorts to ensure high > low order of frequencies
    for key in sorted(sortingdictionary)[::-1]:
        modes.append(sortingdictionary[key])
        freqs.append(key)


    return modes, freqs

    
def printfrequencies(molecule, frequencies, normalmodes):
    bigstring = ""
    for idx,frequency in enumerate(frequencies):
        outstring = ""
        outstring += str(molecule.natom) + "\n"
        
        outstring += "%4.2f" % frequency.real +" cm^-1" + "\n"
        i = 0
        for atom in range(0, molecule.natom):
            outstring += '{: 012.11f}'.format(molecule.geom[atom][0]) + "   "
            outstring += '{: 012.11f}'.format(molecule.geom[atom][1]) + "   "
            outstring += '{: 012.11f}'.format(molecule.geom[atom][2]) + "   "
            outstring += '{: 012.11f}'.format(normalmodes[idx][i])    + "   "
            outstring += '{: 012.11f}'.format(normalmodes[idx][i+1])  + "   "
            outstring += '{: 012.11f}'.format(normalmodes[idx][i+2])  + "   "

            outstring += '\n'
            
            i += 3
            i = i % 8
        
        bigstring += outstring + '\n'
    f = open('example.xyz','w')
    f.write(bigstring)
    f.close()
