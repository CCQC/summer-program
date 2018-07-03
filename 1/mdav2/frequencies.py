import sys
import numpy as np
sys.path.append('../../0/mdav2')
import molecule
import math

ccms = 2.99792458E10 #speed of light in centimeters per second
eh = 4.35974E-18 #J/hartree
ao = 5.29177E-11 #m/br
u = 1.66054E-27 #kg/amu

def readhessian(filename):
    #opens and reads a hessian matrix from the given filename
    with open(filename, 'r') as hessianfile:
        hessianlines = hessianfile.readlines()


    hessianmatrix = np.array([[float(x) for x in hessianlines[line].split()] for line in range(len(hessianlines))])

    return hessianmatrix

def frequencies(molecule,hessian):
    #calculates normal modes and frequencies from a Hessian matrix. Takes a molecule file and a \
    # 3*natom X 3*natom hessian numpy array as input.
    weightedhessian = np.zeros(hessian.shape)
    uwvectors = np.zeros(hessian.shape)
    frequencies = []
    sortedfreqs = []
    normalmodes = []

    #TODO:find a better way to make this matrix
    for A in range(0,molecule.natom * 3):
        for B in range(0, molecule.natom * 3):
            weightedhessian[A][B] = hessian[A][B]/math.sqrt(molecule.masses[math.floor(A/3.)]*molecule.masses[math.floor(B/3.)])

    evals,evecs = np.linalg.eigh(weightedhessian)

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
        if freq.imag:
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
    fout = ""
    for idx, frequency in enumerate(frequencies):
        fout += '{}\n{:4.2f} cm^-1 \n'.format(*[molecule.natom,frequency.real])
        reshape = np.reshape(normalmodes[idx], (3,3))
        for atom, geom, disp in zip(molecule.labels,molecule.geom,reshape):
            lineout = np.concatenate((geom,disp), axis=0).tolist()
            fout += ('{:6}'+'{: 12.11f}   '*6).format(atom,*lineout)+'\n'
        fout += '\n'

    return(fout)
