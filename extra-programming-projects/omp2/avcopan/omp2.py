import psi4.core
import numpy as np
import scipy.linalg as la
from permutation import P

# get global options
MAXITER = psi4.core.get_global_option('MAXITER')
E_CONV  = psi4.core.get_global_option('E_CONVERGENCE')

class SpinOrbitals(object):

  def __init__(self, uhf):
    self.H = uhf.h
    self.G = uhf.g
    self.C = uhf.C
    self.nocc = uhf.nocc            # number of occupied spin-orbitals
    self.norb = uhf.norb            # total number of spin-orbitals = 2 * number of basis functions
    self.nvir = uhf.norb - uhf.nocc # number of virtual spin-orbitals
    self.o = slice(None, self.nocc)
    self.v = slice(self.nocc, None)
    self.mul = lambda T, M : np.tensordot(T, M, axes=(0,0))

  def get_mo_1e_core_hamiltonian(self):
    return self.mul(self.mul(self.H, self.C), self.C)

  def get_mo_2e_repulsion(self):
    return self.mul(self.mul(self.mul(self.mul(self.G, self.C), self.C), self.C), self.C)

  def get_mo_1e_fock_matrix(self):
    o = self.o
    h = self.get_mo_1e_core_hamiltonian()
    g = self.get_mo_2e_repulsion()
    return h + np.einsum('piqi->pq', g[:,o,:,o])

  def rotate(self, U):
    self.C = self.C.dot(U)

class OMP2(object):

  def __init__(self, uhf):
    self.spinorb = SpinOrbitals(uhf)
    self.Vnu = uhf.Vnu
    self.E = 0.0

  def compute_energy(self):
    nocc = self.spinorb.nocc
    nvir = self.spinorb.nvir
    norb = self.spinorb.norb
    t  = np.zeros((nocc, nocc, nvir, nvir))
    o  = self.spinorb.o
    v  = self.spinorb.v
    x  = np.newaxis

    X = np.zeros((norb, norb))
    refG1 = np.zeros((norb, norb))
    corG1 = np.zeros((norb, norb))
    corG2 = np.zeros((norb, norb, norb, norb))
    refG1[o,o] = np.identity(nocc)

    for i in range(MAXITER):
      h = self.spinorb.get_mo_1e_core_hamiltonian()
      g = self.spinorb.get_mo_2e_repulsion()
      f = self.spinorb.get_mo_1e_fock_matrix()
      e = f.diagonal().copy()
      np.fill_diagonal(f, 0.0)

      t = g[o,o,v,v]                                                    \
        + P((2,3)) * np.einsum('ac,ijcb->ijab', f[v,v], t)              \
        - P((0,1)) * np.einsum('ki,kjab->ijab', f[o,o], t)

      t /= (e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v])

      corG1[v,v] = + 1./2 * np.einsum('ijac,ijbc->ab', t, t)
      corG1[o,o] = - 1./2 * np.einsum('jkab,ikab->ij', t, t)
      corG2[o,o,v,v] = t
      corG2[v,v,o,o] = t.T

      G1 = corG1 + refG1
      G2 = corG2                                                        \
         + P((0,1),(2,3)) * np.einsum('pr,qs->pqrs', corG1, refG1)      \
         + P((2,3)) * np.einsum('pr,qs->pqrs', refG1, refG1)

      F = np.einsum('pr,qr->pq', h, G1)                                 \
        + 1./2 * np.einsum('prst,qrst->pq', g, G2)

      X[o,v] = (F - F.T)[o,v] / (e[o,x] - e[x,v])

      U = la.expm(X - X.T)

      self.spinorb.rotate(U)

      E = self.Vnu                                                      \
        + np.einsum('pq,pq', h, G1)                                     \
        + 1./4 * np.einsum('pqrs,pqrs', g, G2)

      dE = E - self.E
      self.E = E

      print         ('@OMP2 {:<3d} {:20.15f} {:20.15f}'  .format(i, E, dE)) # print progress to terminal
      psi4.core.print_out('@OMP2 {:<3d} {:20.15f} {:20.15f}\n'.format(i, E, dE)) # print progress to output
      if(np.fabs(dE) < E_CONV): break  # quit if converged

    return self.E

