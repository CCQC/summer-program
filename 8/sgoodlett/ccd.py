from geom import molecule
from scf import SCF
import numpy as np

class ccd():

    def __init__(self,scf,maxiter=50,conv=6):
        self.scf = scf
        self.E, self.C, self.eps, self.g = self.scf.energy(return_integrals=True)
        self.maxiter = maxiter
        self.conv = 10 ** (-1 * conv)
        self.n_occ = self.scf.n_occ
        self.norbitals = self.scf.norbitals
        self.n_virt = self.norbitals - self.n_occ
        self.t = np.zeros((self.n_occ,self.n_occ,self.n_virt,self.n_virt))
        #self.eps_tensor = np.zeros((self.n_occ,self.n_occ,self.n_virt,self.n_virt))
        #for i in range(self.n_occ):
        #    for j in range(self.n_occ):
        #        for a in range(self.n_virt):
        #            for b in range(self.n_virt):
        #                self.eps_tensor[i,j,a,b] += (self.eps[i] + self.eps[j] - self.eps[a+self.n_occ] - self.eps[b+self.n_occ])**(-1)

    def energy(self, verbose):
        old_energy = 0.0
        round = 0
        self.g = self.g - self.g.transpose((0,3,2,1))
        while True:

            #Update amplitudes (eq. 2)
            for i in range(self.n_occ):
                for j in range(self.n_occ):
                    for a in range(self.n_occ, self.norbitals):
                        for b in range(self.n_occ, self.norbitals):
                            ap = a - self.n_occ
                            bp = b - self.n_occ
                            Matthew = 0
                            Mark = 0
                            Luke = 0
                            John = 0
                            #All the stuff inside the parentheses
                            for c in range(self.n_occ,self.norbitals):
                                for d in range(self.n_occ,self.norbitals):
                                    Matthew += 0.5 * self.g[a,b,c,d] * self.t[i,j,c-self.n_occ,d-self.n_occ]

                            for k in range(self.n_occ):
                                for l in range(self.n_occ):
                                    Mark += 0.5 * self.g[k,l,i,j] * self.t[k,l,ap,bp]

                            for k in range(self.n_occ):
                                for c in range(self.n_occ,self.norbitals):
                                    cp = c - self.n_occ
                                    Luke += (self.g[a,k,i,c] * self.t[j,k,bp,cp]) - (self.g[b,k,i,c] * self.t[j,k,ap,cp]) + (self.g[b,k,j,c] * self.t[i,k,ap,cp]) - (self.g[a,k,j,c] * self.t[i,k,bp,cp])

                            for k in range(self.n_occ):
                                for l in range(self.n_occ):
                                    for c in range(self.n_virt):
                                        for d in range(self.n_virt):
                                            cp = c + self.n_occ
                                            dp = d + self.n_occ
                                            #Like the Gospel of John, this is going to be very long and involved
                                            Verbum = -0.5 * ((self.t[i,j,ap,c] * self.t[k,l,bp,d]) - (self.t[i,j,bp,c] * self.t[k,l,ap,d]))
                                            caro = -0.5 * ((self.t[i,k,ap,bp] * self.t[j,l,c,d]) - (self.t[j,k,ap,bp] * self.t[i,l,c,d]))
                                            factum = 0.25 * self.t[i,j,c,d] * self.t[k,l,ap,bp]
                                            est = (self.t[i,k,ap,c] * self.t[j,l,bp,d]) - (self.t[j,k,ap,c] * self.t[i,l,bp,d])
                                            John += self.g[k,l,cp,dp] * (Verbum + caro + factum + est)

                            eps = (self.eps[i] + self.eps[j] - self.eps[a] - self.eps[b])**(-1)
                            self.t[i,j,ap,bp] = eps * (self.g[a,b,i,j] + Matthew + Mark + Luke + John)

            #Evaluate correlation energy (eq. 1)
            energy = 0
            for k in range(self.n_occ):
                for l in range(self.n_occ):
                    for c in range(self.n_virt):
                        for d in range(self.n_virt):
                            energy += self.g[k,l,c+self.n_occ,d+self.n_occ]*self.t[k,l,c,d]
            energy *= 0.25

            if verbose:
                print("CCD Correlation Energy: {} Hartrees".format(energy))
            #Convergence check
            if abs(energy - old_energy) < self.conv:
                return energy
            if round >= self.maxiter:
                print("Surpassed max CC iteration")
                break

            old_energy = energy
            round += 1

if __name__ == '__main__':
    mol = molecule(units='Bohr')
    mol.to_angstrom()
    scf = SCF(mol,basis='STO-3G',maxiter=150)
    ccd = ccd(scf)
    ccd.energy(verbose=True)