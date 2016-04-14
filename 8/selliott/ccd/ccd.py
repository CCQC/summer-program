import psi4
import numpy as np
from perm import P

class CCD(object):

  def __init__(self, uhf):
    
    self.norb = uhf.norb            #number of orbitals
    self.nocc = uhf.nocc            #number of occupied orbitals 
    self.nvir = uhf.norb-uhf.nocc   #number of virtual orbitals
    
    #Factorizing the AO to MO transformation (following Sherill's DF ERI notes)
    G = uhf.g #antisymmetrized AO-basis TEI <mu nu||rho sigma>
    C = uhf.C #MO coefficients
    G = np.einsum('Pu,Pvrs->uvrs',C,G) 
    G = np.einsum('Qv,PQrs->Pvrs',C,G)
    G = np.einsum('Rr,PQRs->PQrs',C,G)
    G = np.einsum('Ss,PQRS->PQRs',C,G)
 
    self.g=G
    self.e = uhf.e #orbital energies
    self.tstart = np.zeros((self.nocc,self.nocc,self.nvir,self.nvir)) #initial amplitude
    self.E = 0.0 #initial energy
 
  def update_energy(self,t):
    """
    Finds the CCD correlation energy
    E = 1/4 Sum(klcd) <kl || cd >>t_kl^cd
    """
    o = slice(None, self.nocc)
    v = slice(self.nocc, None)
    Ec = 1/4. * np.sum(self.g[o,o,v,v]*t)
    
    return Ec

  def update_amplitudes(self,t):
    """
    updates the amplitues
    """
    o = slice(None, self.nocc)
    v = slice(self.nocc, None)
    x = np.newaxis
    g = self.g

    HN = g[o,o,v,v]

    HNT2 = 1/2.* np.einsum('abcd,ijcd->ijab', g[v,v,v,v], t) 
    HNT2 += 1/2.* np.einsum('klij,klab->ijab',g[o,o,o,o], t)
    HNT2 += P(0,1) *( P(2,3) * np.einsum('akic,jkbc->ijab',g[v,o,o,v],t) )

    HNT2T2 =  -1/2. * (P(2,3) * np.einsum('klcd,ijac,klbd->ijab',g[o,o,v,v],t,t))
    HNT2T2 += -1/2. * (P(0,1) * np.einsum('klcd,ikab,jlcd->ijab',g[o,o,v,v],t,t))
    HNT2T2 += 1/4.  * np.einsum('klcd,ijcd,klab->ijab',g[o,o,v,v],t,t)
    HNT2T2 += P(0,1)* np.einsum('klcd,ikac,jlbd->ijab',g[o,o,v,v],t,t)

    eps = 1.0 / (self.e[o,x,x,x] + self.e[x,o,x,x] - self.e[x,x,v,x] - self.e[x,x,x,v])
    t = (HN + HNT2 +  HNT2T2) * eps
    
    return t

  def corr_energy(self):

    told = self.tstart
    Eold = self.E
    for i in range(50):
      
      tnew = self.update_amplitudes(told)
      Enew = self.update_energy(tnew)
      print( "CCD " + str(i)+" \t"+str(Enew)+" \t"+" \t"+str(Enew-Eold) )
      
      #if(Enew - Eold) < Econv:
      if(abs(Enew - Eold) < 0.000000000001):
        break

      Eold = Enew
      told = tnew
       
    return Enew
       
        
