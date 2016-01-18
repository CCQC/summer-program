import psi4
import numpy as np
from permutation import P

# get global options
MAXITER = psi4.get_global_option('MAXITER')
E_CONV  = psi4.get_global_option('E_CONVERGENCE')

class CCD(object):

  def __init__(self, uhf):
    G = uhf.g # antisymmetrized AO-basis TEIs, <mu nu || rh si>
    C = uhf.C # MO coefficients
    self.e    = uhf.e               # orbital energies
    self.nocc = uhf.nocc            # number of occupied spin-orbitals
    self.norb = uhf.norb            # total number of spin-orbitals = 2 * number of basis functions
    self.nvir = uhf.norb - uhf.nocc # number of virtual spin-orbitals

    # transform integrals to MO basis
    mul = lambda T, M : np.tensordot(T, M, axes=(0,0))
    self.g = mul(mul(mul(mul(G, C), C), C), C)

    self.E = 0.0

  def compute_energy(self):
    nocc, nvir, e, g = self.nocc, self.nvir, self.e, self.g
    t  = np.zeros((nocc, nocc, nvir, nvir))
    o  = slice(None, nocc)
    v  = slice(nocc, None)
    x = np.newaxis
    Ep = 1./(e[o,x,x,x] + e[x,o,x,x] - e[x,x,v,x] - e[x,x,x,v])

    for i in range(MAXITER):
      # update T2 amplitudes
      t  = g[o,o,v,v]                                                                  \
         + 1./2                  * np.einsum("abcd,ijcd->ijab", g[v,v,v,v], t)         \
         + 1./2                  * np.einsum("klij,klab->ijab", g[o,o,o,o], t)         \
         + 1.   * P((0,1),(2,3)) * np.einsum("akic,jkbc->ijab", g[v,o,o,v], t)         \
         - 1./2 * P((2,3))       * np.einsum("klcd,ijac,klbd->ijab", g[o,o,v,v], t, t) \
         - 1./2 * P((0,1))       * np.einsum("klcd,ikab,jlcd->ijab", g[o,o,v,v], t, t) \
         + 1./4                  * np.einsum("klcd,ijcd,klab->ijab", g[o,o,v,v], t, t) \
         + 1.   * P((0,1))       * np.einsum("klcd,ikac,jlbd->ijab", g[o,o,v,v], t, t)
      t *= Ep
      # evaluate energy
      E  = 1./4 * np.sum(g[o,o,v,v] * t)
      dE = E - self.E
      self.E, self.t = E, t
      print         ('@CCD {:<3d} {:20.15f} {:20.15f}'  .format(i, E, dE)) # print progress to terminal
      psi4.print_out('@CCD {:<3d} {:20.15f} {:20.15f}\n'.format(i, E, dE)) # print progress to output
      if(np.fabs(dE) < E_CONV): break  # quit if converged

    return self.E

