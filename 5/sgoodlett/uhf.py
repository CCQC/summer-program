from scf import SCF
import numpy as np
from geom import molecule
from scipy.linalg import inv, sqrtm, eigh, block_diag

class UHF(SCF):

    def __init__(self,mol,basis,convergence=6,maxiter=50):

        SCF.__init__(self,mol,basis,convergence=convergence,maxiter=maxiter)
        n2 = self.norbitals
        self.norbitals = 2 * self.norbitals
        self.n_occ = self.n_e
        self.S = block_diag(self.S,self.S)
        self.T = block_diag(self.T,self.T)
        self.V = block_diag(self.V,self.V)
        self.g.transpose((0,2,1,3))
        self.g_ = np.zeros((self.norbitals,self.norbitals,self.norbitals,self.norbitals))
        G = np.zeros((n2,n2,self.norbitals,self.norbitals))
        for mu in range(n2):
            for nu in range(n2):
                G[mu,nu] += block_diag(self.g_chem[mu,nu,:,:],self.g_chem[mu,nu,:,:])
        self.g_[0:n2,0:n2,:,:] += G
        self.g_[n2:,n2:,:,:] += G

        self.g = self.g_
        self.g_chem = self.g
        self.gt = self.g.transpose((0,3,2,1))
        self.X = np.matrix(inv(sqrtm(self.S)))
        self.D = np.zeros((self.norbitals,self.norbitals))

    def scf_cycle(self):
        #Build Fock matrix (eq 4)
        v = np.einsum('mnrs,sr->mn',self.g - self.gt, self.D)
        h = self.T + self.V
        f = h + v

        #Transform f -> ~f to orthogonalized AO basis (eq 6)
        f_ortho = np.matmul(self.X,np.matmul(f,self.X))

        #Diagonalize ~f, yielding orbital energies \epsilon_p and MO coeffecients ~C_{\mu p} (eq. 6)
        eps, C_ortho = eigh(f_ortho)

        #Backtransform ~C -> C to original AO basis (eq. 6)
        C = np.matmul(self.X,np.matrix(C_ortho))

        #Build Density matrix D (eq. 3)
        self.D = np.einsum('ij,kj->ik',C[:,:self.n_occ],np.conj(C[:,:self.n_occ]))

        #Evaluate Energy (eq. 5)
        hv = h + (0.5 * v)
        energy = np.einsum('ij,ji->',hv,self.D) + self.E_nuc

        return energy, C, eps


if __name__ == '__main__':
    mol = molecule(units='Bohr')
    mol.to_angstrom()
    E = UHF(mol,'STO-3G')
    e = E.energy(verbose=True)

