import psi4
import numpy as np

class CCD:

  def __init__(self, uhf):
    
    self.norb = uhf.norb            #number of orbitals
    self.nocc = uhf.nocc            #number of occupied orbitals 
    self.nvir = uhf.norb-uhf.nocc   #number of virtual orbitals
    
    #Factorizing the AO to MO transformation (following Sherill's DF ERI notes)
    G = uhf.g #antisymmetrized AO-basis TEI <mu nu||rho sigma>
    C = uhf.C #MO coefficients
    G = np.einsum('Pu','Pvrs->uvrs',C,G) 
    G = np.einsum('Qv','PQrs->Pvrs',C,G)
    G = np.einsum('Rr','PQRs->PQrs',C,G)
    G = np.einsum('Ss','PQRS->PQRs',C,G)

    self.g=G
    self.e = uhf.e #orbital energies
    self.t = np.zeros((self.nocc,self.nocc,self.nvir,self.nvir)) #initial amplitude
    self.E = 0.0 #initial energy

    def energy(self):
    """
    Finds the CCD correlation energy
    E = 1/4 Sum(klcd) <kl || cd >>t_kl^cd
    """
    Ec = 1/4 * np.sum(self.g[self.nocc,self.nocc,self.vir,self.vir]*self.t)
    return Ec

    def amplitudes(self):

