import psi4
import numpy as np
from perm import P

class CCSD(object):

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
    self.dstart = np.zeros((self.nocc,self.nocc,self.nvir,self.nvir)) #initial doubles amplitude
    self.sstart = np.zeros((self.nocc,self.nvir)) #initial singles amplitude
    self.E = 0.0 #initial energy
    
    
    self.o = o = slice(None, self.nocc)
    self.v = v = slice(self.nocc, None)
    self.x = x = np.newaxis
    self.Eps = 1.0 / (self.e[o,x,x,x] + self.e[x,o,x,x] - self.e[x,x,v,x] - self.e[x,x,x,v])
    self.eps = 1.0 / (self.e[o,x] - self.e[x,v])
    self.dstart = G[o,o,v,v]*self.Eps
   
  def update_energy(self,S,D):
    """
    Finds the CCSD correlation energy
    """
    o,v =self.o,self.v 
    Ec = 1/4. * np.einsum('ijab,ijab',self.g[o,o,v,v],D) +  1/2. * np.einsum('ijab,ia,jb',self.g[o,o,v,v],S,S) 
    return Ec

  
  def update_singles(self,S,D):
    """
    updates the singles amplitues
    """
    g,o,v,x,eps = self.g,self.o,self.v,self.x,self.eps

    s = eps * (
       1/2.* np.einsum('akcd,ikcd->ia',g[v,o,v,v],D)        \
      -1/2.* np.einsum('klic,klac->ia',g[o,o,o,v],D)          \
      +      np.einsum('akic,kc->ia',g[v,o,o,v],S)            \
      -1/2.* np.einsum('klcd,ic,klad->ia',g[o,o,v,v],S,D)     \
      -1/2.* np.einsum('klcd,ka,ilcd->ia',g[o,o,v,v],S,D)     \
      +      np.einsum('klcd,kc,lida->ia',g[o,o,v,v],S,D)     \
      +      np.einsum('akcd,ic,kd->ia',g[v,o,v,v],S,S)       \
      -      np.einsum('klic,ka,lc->ia',g[o,o,o,v],S,S)       \
      -      np.einsum('klcd,ic,ka,ld->ia',g[o,o,v,v],S,S,S) )

    return s 


  def update_doubles(self,S,D):
    """
    updates the amplitues
    """
    g,o,v,x,Eps = self.g,self.o,self.v,self.x,self.Eps

    #singles amplitude contribution
    d = Eps * (g[o,o,v,v]
         + 1/2.*                     np.einsum('abcd,ijcd->ijab', g[v,v,v,v], D)        \
         + 1/2.*                     np.einsum('klij,klab->ijab',g[o,o,o,o], D)         \
         +        P(0,1) *( P(2,3) * np.einsum('akic,jkbc->ijab',g[v,o,o,v],D) )        \
         - 1/2.* (P(2,3) *           np.einsum('klcd,ijac,klbd->ijab',g[o,o,v,v],D,D))  \
         - 1/2.* (P(0,1) *           np.einsum('klcd,ikab,jlcd->ijab',g[o,o,v,v],D,D))  \
         + 1/4.*                     np.einsum('klcd,ijcd,klab->ijab',g[o,o,v,v],D,D)   \
         +        P(0,1) *           np.einsum('klcd,ikac,jlbd->ijab',g[o,o,v,v],D,D)  \
         +        P(0,1) *           np.einsum('abcj,ic->ijab',g[v,v,v,o],S) \
         -       (P(2,3) *           np.einsum('kbij,ka->ijab',g[o,v,o,o],S))\
         +        P(0,1) * (P(2,3) * np.einsum('akcd,ic,kjdb->ijab',g[v,o,v,v],S,D) )\
         -       (P(0,1) * (P(2,3) * np.einsum('klic,ka,ljcb->ijab',g[o,o,o,v],S,D) ) ) \
         - 1/2.* (P(2,3) *           np.einsum('kbcd,ka,ijcd->ijab',g[o,v,v,v],S,D))              \
         + 1/2.* (P(0,1) *           np.einsum('klcj,ic,klab->ijab',g[o,o,v,o],S,D))             \
         +        P(2,3) *           np.einsum('kacd,kc,ijdb->ijab',g[o,v,v,v],S,D)                    \
         -       (P(0,1) *           np.einsum('klci,kc,ljab->ijab',g[o,o,v,o],S,D) )                  \
         +                           np.einsum('abcd,ic,jd->ijab',g[v,v,v,v],S,S)                                  \
         +                           np.einsum('klij,ka,lb->ijab',g[o,o,o,o],S,S)                                  \
         -       (P(0,1) * (P(2,3) * np.einsum('kbcj,ic,ka->ijab',g[o,v,v,o],S,S)) )         \
         + 1/2.*                     np.einsum('klcd,ic,jd,klab->ijab',g[o,o,v,v],S,S,D)                    \
         + 1/2.*                     np.einsum('klcd,ka,lb,ijcd->ijab',g[o,o,v,v],S,S,D)                    \
         -       (P(0,1) * (P(2,3) * np.einsum('klcd,ic,ka,ljdb->ijab',g[o,o,v,v],S,S,D)) ) \
         -       (P(0,1) *           np.einsum('klcd,kc,id,ljab->ijab',g[o,o,v,v],S,S,D))              \
         -       (P(2,3) *           np.einsum('klcd,kc,la,ijdb->ijab',g[o,o,v,v],S,S,D))              \
         +        P(2,3) *           np.einsum('kbcd,ic,ka,jd->ijab',g[o,v,v,v],S,S,S)                    \
         +        P(0,1) *           np.einsum('klcj,ic,ka,lb->ijab',g[o,o,v,o],S,S,S) \
         +                           np.einsum('klcd,ic,jd,ka,lb->ijab',g[o,o,v,v],S,S,S,S) )
    return d

  def corr_energy(self):
    """
    CCSD iterations
    """
    Sold = self.sstart
    Dold = self.dstart
    Eold = self.E
    for i in range(50):
      Snew = self.update_singles(Sold,Dold)
      Dnew = self.update_doubles(Snew,Dold)
      Enew = self.update_energy(Snew,Dnew)
      print( "CCSD " + str(i)+" \t"+str(Enew)+" \t"+" \t"+str(Enew-Eold) )
      
      if(abs(Enew - Eold) < 0.000000000001):
        break

      Sold = Snew
      Dold = Dnew
      Eold = Enew
       
    return Enew
       
        
