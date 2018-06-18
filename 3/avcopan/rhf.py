import psi4.core
import numpy        as np
import scipy.linalg as la
from tensor import tensor # subclass of ndarray with operator% overloaded for matrix multiplication

# get global options
MAXITER = psi4.core.get_global_option('MAXITER')
E_CONV  = psi4.core.get_global_option('E_CONVERGENCE')

class RHF(object):

  def __init__(self, mol, mints):
    self.mol   = mol
    self.mints = mints
    self.E   = 0.0
    self.Vnu = mol.nuclear_repulsion_energy()

    # one-electron integrals
    S = np.array(mints.ao_overlap()).view(tensor)
    T = np.array(mints.ao_kinetic()).view(tensor)
    V = np.array(mints.ao_potential()).view(tensor)
    self.h = T + V # core Hamiltonian
    self.X = la.inv(la.sqrtm(S)).view(tensor) # orthogonalizer
    # two-electron integrlas, transposed to physicist's notation
    self.g = np.array(mints.ao_eri()).transpose((0,2,1,3))

    # get totla number of spatial orbitals
    self.norb = mints.basisset().nbf()
    # get number of occupied spatial orbitals
    nelec = sum(mol.Z(A) for A in range(mol.natom())) - mol.molecular_charge() # num e- = sum(Z_A) - mol. charge
    self.nocc = int(nelec/2)

    # if it isn't a closed-shell molecule, complain
    if mol.multiplicity() != 1 or nelec % 2 != 0:
      raise Exception("RHF code only takes closed-shell molecules")

  def compute_energy(self):
    norb, nocc, Vnu, X, h, g = self.norb, self.nocc, self.Vnu, self.X, self.h, self.g

    # start with D = 0 for core guess
    D = np.zeros((norb, norb)).view(tensor)

    for i in range(MAXITER):
      v  =        np.tensordot(g, D, axes = [(1,3), (1,0)]) # v_mu,nu  =        sum_rh,si <mu rh|nu si> D_si,rh
      v -= 1./2 * np.tensordot(g, D, axes = [(1,2), (1,0)]) # v_mu,nu -= 1./2 * sum_rh,si <mu rh|si nu> D_si,rh
      f  = h + v                       # build Fock matrix

      tf    = X % f % X                # transform to orthogonalized AO basis
      e, tC = la.eigh(tf)              # diagonalize
      C     = X % tC                   # backtransform to original AO basis
      oC    = C[:,:nocc]               # chop out occupied block of MO coefficient matrix
      D     = 2 * oC % oC.T            # update density matrix

      E = np.trace((h + v/2) % D) + Vnu  # compute energy by tracing with density matrix

      dE = E - self.E                  # get change in energy
      self.E, self.C, self.e = E, C, e # save these for later
      psi4.core.print_out('@RHF {:<3d} {:20.15f} {:20.15f}\n'.format(i, E, dE)) # print progress
      if(np.fabs(dE) < E_CONV): break  # quit if converged

    return self.E
