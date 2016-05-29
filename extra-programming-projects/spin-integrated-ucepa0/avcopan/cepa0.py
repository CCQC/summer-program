import psi4
import numpy as np
from permutation_operator import P

# get global options
MAXITER = psi4.get_global_option('MAXITER')
E_CONV  = psi4.get_global_option('E_CONVERGENCE')

class CEPA0(object):

  def __init__(self, uhf):
    G  = uhf.g
    aC, bC = uhf.aC, uhf.bC
    ae, be = uhf.ae, uhf.be
    nbf   = uhf.nbf
    naocc = uhf.naocc
    nbocc = uhf.nbocc
    navir = nbf - naocc
    nbvir = nbf - nbocc

    # built guess (empty) t-amplitudes
    aaT = np.zeros((naocc, naocc, navir, navir))
    abT = np.zeros((naocc, nbocc, navir, nbvir))
    bbT = np.zeros((nbocc, nbocc, nbvir, nbvir))

    # built resolvent denominators
    ao, bo = slice(None, naocc), slice(None, nbocc)
    av, bv = slice(naocc, None), slice(nbocc, None)
    x = np.newaxis
    aaEp = 1./(ae[ao,x,x,x] + ae[x,ao,x,x] - ae[x,x,av,x] - ae[x,x,x,av])
    abEp = 1./(ae[ao,x,x,x] + be[x,bo,x,x] - ae[x,x,av,x] - be[x,x,x,bv])
    bbEp = 1./(be[bo,x,x,x] + be[x,bo,x,x] - be[x,x,bv,x] - be[x,x,x,bv])

    # transform integrals to MO basis
    mul = lambda T, M : np.tensordot(T, M, axes=(0,0))
    aaG = mul(mul(mul(mul(G, aC), aC), aC), aC)
    abG = mul(mul(mul(mul(G, aC), bC), aC), bC)
    bbG = mul(mul(mul(mul(G, bC), bC), bC), bC)
    aaG = aaG - aaG.transpose((0,1,3,2)) # antisymmetrize diagonal blocks
    bbG = bbG - bbG.transpose((0,1,3,2))

    # save what we need
    self.ao, self.bo, self.av, self.bv = ao, bo, av, bv
    self.aaT , self.abT , self.bbT  = aaT , abT , bbT
    self.aaG , self.abG , self.bbG  = aaG , abG , bbG
    self.aaEp, self.abEp, self.bbEp = aaEp, abEp, bbEp
    self.E = 0.0

  def compute_energy(self):
    ao, bo, av, bv = self.ao, self.bo, self.av, self.bv
    aaT , abT , bbT  = np.zeros(self.aaT.shape), np.zeros(self.abT.shape), np.zeros(self.bbT.shape)
    aaG , abG , bbG  = self.aaG , self.abG , self.bbG 
    aaEp, abEp, bbEp = self.aaEp, self.abEp, self.bbEp

    for i in range(MAXITER):
      self.aaT = aaEp * (  aaG[ao,ao,av,av]                                                   \
                         + 1./2         * np.einsum("abcd,ijcd->ijab", aaG[av,av,av,av], aaT) \
                         + 1./2         * np.einsum("klij,klab->ijab", aaG[ao,ao,ao,ao], aaT) \
                         + P("0/1|2/3") * np.einsum("akic,jkbc->ijab", aaG[av,ao,ao,av], aaT) \
                         + P("0/1|2/3") * np.einsum("aKiC,jKbC->ijab", abG[av,bo,ao,bv], abT) )

      self.abT = abEp * (  abG[ao,bo,av,bv]                                    \
                         + np.einsum("aBcD,iJcD->iJaB", abG[av,bv,av,bv], abT) \
                         + np.einsum("kLiJ,kLaB->iJaB", abG[ao,bo,ao,bo], abT) \
                         + np.einsum("akic,kJcB->iJaB", aaG[av,ao,ao,av], abT) \
                         + np.einsum("aKiC,JKBC->iJaB", abG[av,bo,ao,bv], bbT) \
                         - np.einsum("aKcJ,iKcB->iJaB", abG[av,bo,av,bo], abT) \
                         - np.einsum("kBiC,kJaC->iJaB", abG[ao,bv,ao,bv], abT) \
                         + np.einsum("kBcJ,ikac->iJaB", abG[ao,bv,av,bo], aaT) \
                         + np.einsum("BKJC,iKaC->iJaB", bbG[bv,bo,bo,bv], abT) )

      self.bbT = bbEp * (  bbG[bo,bo,bv,bv]                                                   \
                         + 1./2         * np.einsum("ABCD,IJCD->IJAB", bbG[bv,bv,bv,bv], bbT) \
                         + 1./2         * np.einsum("KLIJ,KLAB->IJAB", bbG[bo,bo,bo,bo], bbT) \
                         + P("0/1|2/3") * np.einsum("AKIC,JKBC->IJAB", bbG[bv,bo,bo,bv], bbT) \
                         + P("0/1|2/3") * np.einsum("kAcI,kJcB->IJAB", abG[ao,bv,av,bo], abT) )

      aaT, abT, bbT = self.aaT, self.abT, self.bbT

      # evaluate the correlation energy
      E = 1./4 * np.sum(aaG[ao,ao,av,av] * aaT) \
        + 1.   * np.sum(abG[ao,bo,av,bv] * abT) \
        + 1./4 * np.sum(bbG[bo,bo,bv,bv] * bbT)

      dE = E - self.E
      self.E = E
      print         ('@CEPA0 {:<3d} {:20.15f} {:20.15f}'  .format(i, E, dE)) # print progress to terminal
      psi4.print_out('@CEPA0 {:<3d} {:20.15f} {:20.15f}\n'.format(i, E, dE)) # print progress to output
      if(np.fabs(dE) < E_CONV): break  # quit if converged

    return E
