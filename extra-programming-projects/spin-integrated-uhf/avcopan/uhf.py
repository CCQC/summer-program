import psi4
import numpy as np
import scipy.linalg as la

class UHF(object):

  def __init__(self, mol, mints):
    nbf = mints.basisset().nbf()
    S = np.matrix(mints.ao_overlap())
    T = np.matrix(mints.ao_kinetic())
    V = np.matrix(mints.ao_potential())
    G = np.array(mints.ao_eri())
    # object attributes
    self.h  = T + V
    self.g  = G.transpose(0,2,1,3) # un-antisymmetrized phys. notation integrals
    self.aD = np.zeros((nbf, nbf))
    self.bD = np.zeros((nbf, nbf))
    self.X  = np.matrix(la.inv(la.sqrtm(S)))
    self.vnu = mol.nuclear_repulsion_energy()
    # work out the number of alpha and beta electrons, assuming high-spin open-shell
    nelec = int( sum(mol.Z(A) for A in range(mol.natom())) - mol.molecular_charge() ) # num e- = sum(Z_A) - mol. charge
    nsocc = mol.multiplicity() - 1 # number of unpaired eletrons = 2 * S
    ndocc = (nelec - nsocc) / 2    # number of doubly occupied orbitals = number of electrons - number unpaired
    self.naocc = ndocc + nsocc     # number of alpha electrons
    self.nbocc = ndocc             # number of beta electrons
    self.nbf   = nbf

  def compute_energy(self):
    h, g, aD, bD, X, vnu, naocc, nbocc = self.h, self.g, self.aD, self.bD, self.X, self.vnu, self.naocc, self.nbocc
    self.E = 0.0

    for i in range(psi4.get_option('SCF', 'MAXITER')):
      aV = np.einsum('mrns,sr', g, aD + bD) - np.einsum('mrsn,sr', g, aD)
      bV = np.einsum('mrns,sr', g, aD + bD) - np.einsum('mrsn,sr', g, bD)
      aF = h + aV
      bF = h + bV
      taF = X * aF * X
      tbF = X * bF * X
      ae, taC = la.eigh(taF)
      be, tbC = la.eigh(tbF)
      aC  = X * taC
      bC  = X * tbC
      oaC = aC[:,:naocc]
      obC = bC[:,:nbocc]
      aD  = oaC * oaC.T
      bD  = obC * obC.T
      E = np.trace((h + aV/2)*aD) + np.trace((h + bV/2)*bD) + vnu
      dE = E - self.E
      self.E, self.aC, self.bC, self.ae, self.be = E, aC, bC, ae, be
      print('UHF {:-3d} {:20.15f} {:20.15f}'.format(i, E, dE)) # print progress
      if(np.fabs(dE) < psi4.get_global_option('E_CONVERGENCE')): break # quit if converged

