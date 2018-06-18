# An internal coordinate vibrational analaysis class, which takes G and F matrices
# as its constructor argument as well as an optional "blockdims" argument indicating
# the symmetry-blocking of F and G.  For example, blockdims=(4, 2, 3, 3,) indicates
# that there are 4, 2, 3, and 3 symmetrized internal coordinates that transform as the
# first, second, third, and fourth irreps, respectively.  The actual work of the GF
# matrix method happens in the function self.compute_frequencies(), and there are
# two functions for printing the results: self.print_cm_frequencies() prints the
# frequencies (blocked by irrep) in inverse centimeters, and the function
# self.print_total_energy_distribution() prints the contribution of each internal
# coordinate to the normal modes (following the same ordering as the frequencies).
import numpy as np
import scipy.linalg as la
from molecule import Molecule

class IntCoVibAnalysis(object):

  def __init__(self, G, F, blockdims=None):
    self.nintcos, = G.diagonal().shape
    if blockdims is None: blockdims = (self.nintcos,)
    self.blockdims = blockdims # (dim1, dim2, ...) for when G and F are blocked by symmetry
    self.nirreps = len(blockdims)
    self.G = G
    self.F = F
    self.X = la.sqrtm(G) # orthogonalizer
    self.L  = np.zeros((self.nintcos, self.nintcos))
    self.gf = np.zeros((self.nintcos,))
    self.slices  = self.build_block_slice_vector()
    self.vlabels = ["v{:d}".format(n+1) for n in range(self.nintcos)]
    self.slabels = ["s{:d}".format(k+1) for k in range(self.nintcos)]
    self.compute_frequencies()

  def build_block_slice_vector(self):
    starts = [sum(self.blockdims[:i]) for i in range(self.nirreps)]
    ends   = starts[1:] + [None]
    return [slice(start, end) for start, end in zip(starts, ends)]

  def compute_frequencies(self):
    for slc in self.slices:
      X, F = self.X[slc, slc], self.F[slc, slc] # diagonalize by block
      tF = np.dot(X, np.dot(F, X))
      gf, tL = la.eigh(tF)
      L = np.dot(X, tL)
      self.gf[slc], self.L[slc,slc] = gf, L
    aJ2J, amu2kg, ang2m, c = 1e-18, 1.6605389e-27, 1e-10, 29979245800.0 # speed of light in cm/s
    self.frequencies = np.sqrt(self.gf) * np.sqrt(aJ2J/(amu2kg*ang2m*ang2m))/(c*2*np.pi)
    

  def print_cm_frequencies(self):
    for irrep, slc in enumerate(self.slices):
      vlabels = self.vlabels[slc]
      frequencies = self.frequencies[slc][::-1]
      print("\nFrequencies for Irrep Block {:d}:".format(irrep))
      for label, freq in zip(vlabels, frequencies):
        print("{:<3s} {: >7.2f}  cm^-1".format(label, freq))

  def print_total_energy_distribution(self):
    for irrep, slc, in enumerate(self.slices):
      L   = self.L[slc, slc]
      iL  = la.inv(L)
      TED = (L * (iL).T)[:,::-1] # reverse column order for consistent labeling
      slabels = self.slabels[slc]
      vlabels = self.vlabels[slc]
      print("\nTotal Energy Distribution for Irrep Block {:d}:".format(irrep))
      print("{:<3s}".format('') + ''.join(" {:>6s}".format(slabel) for slabel in slabels))
      for n, vlabel in enumerate(vlabels):
        print("{:<3s}".format(vlabel) + ''.join(" {:>6.3f}".format(elem) for elem in TED[:,n]))

