import sys
import numpy as np
sys.path.append('../../extra-files')
sys.path.append('../../0/mdav2')
import molecule
import math

ccms = 3E10 #speed of light in centimeters per second
eh = 4.36E-18 #J/hartree
ao = 5.29E-11 #m/br
u = 1.66E-27 #kg/amu

def readhessian(filename):
    #opens and reads a hessian matrix from the given filename
    with open(filename, 'r') as hessianfile:
        hessianlines = hessianfile.readlines()


    hessianmatrix = np.zeros(shape=(len(hessianlines),len(hessianlines)))
    for line in range(0,len(hessianlines)):
        hessianmatrix[line] = hessianlines[line].split()

    return hessianmatrix

def frequencies(molecule,hessian):
    #calculates normal modes and frequencies from a Hessian matrix. Takes a molecule file and a \
    # 3*natom X 3*natom hessian numpy array as input.
    weightedhessian = np.zeros(hessian.shape)
    uwvectors = np.zeros(weightedhessian.shape)
    frequencies = []
    sortedfreqs = []
    normalmodes = []

    #TODO:find a better way to make this matrix
    for A in range(0,molecule.natom * 3):
        for B in range(0, molecule.natom * 3):
            weightedhessian[A][B] = hessian[A][B]/math.sqrt(molecule.masses[math.floor(A/3.)]*molecule.masses[math.floor(B/3.)])

    evals,evecs = np.linalg.eig(weightedhessian)

    #TODO: work on these lines to a more readable state
    for vector in range(0,molecule.natom * 3):
        uwvectors[vector] = evecs[vector]/math.sqrt(molecule.masses[math.floor(vector/3.)])

    #transpose to make slicing easier
    uwvectors = uwvectors.T

    #convert from atomic units to SI
    evals = evals * (eh/(u * ao **2))

    for idx,frequency in enumerate(evals):
        #TODO: find a neater way to deal with negative force constants/imaginary frequencies
        #convert from force constants to wavenumbers
        freq = ( (1./(math.pi * 2) * np.sqrt(frequency + 0J))/ccms) #+0J to force complex evaluation
        #turns imaginary frequencies to negative
        if freq.imag != 0:
            freq = -1*freq.imag
        frequencies.append(freq)

    #make a dictionary of frequencies and modes to sort frequencies high>low
    sortingdictionary = dict(zip(frequencies,uwvectors))

    #sorts to ensure high > low order of frequencies
    for key in sorted(sortingdictionary)[::-1]:
        normalmodes.append(sortingdictionary[key])
        sortedfreqs.append(key)


    return normalmodes, sortedfreqs


def printfrequencies(molecule, frequencies, normalmodes):
    #TODO: use cleaner formatting tools here
    finaloutput = ""
    for idx,frequency in enumerate(frequencies):
        tempout = ""
        tempout += str(molecule.natom) + "\n"

        tempout += "%4.2f" % frequency.real +" cm^-1" + "\n"
        i = 0
        for atom in range(0, molecule.natom):
            tempout += '{: 012.11f}'.format(molecule.geom[atom][0]) + "   "
            tempout += '{: 012.11f}'.format(molecule.geom[atom][1]) + "   "
            tempout += '{: 012.11f}'.format(molecule.geom[atom][2]) + "   "
            tempout += '{: 012.11f}'.format(normalmodes[idx][i])    + "   "
            tempout += '{: 012.11f}'.format(normalmodes[idx][i+1])  + "   "
            tempout += '{: 012.11f}'.format(normalmodes[idx][i+2])  + "   "

            tempout += '\n'

            i += 3
            i = i % 8

        finaloutput += tempout + '\n'
    with open('example.xyz','w') as f:
        f.write(finaloutput)
