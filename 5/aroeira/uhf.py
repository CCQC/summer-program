import psi4
import numpy as np
import scipy.linalg as la

class UHF:

    def __init__(self, mol, mints):

        self.mol = mol
        self.mints = mints

        self.V_nuc = mol.nuclear_repulsion_energy()
        self.T = np.matrix(mints.ao_kinetic())
        self.S = np.matrix(mints.ao_overlap())
        self.V = np.matrix(mints.ao_potential())

        self.g = np.array(mints.ao_eri()).swapaxes(1,2)

        self.nelec = -mol.molecular_charge()
        for A in range(mol.natom()):
            self.nelec += int(mol.Z(A))
        self.ndocc = int(self.nelec / 2)
        self.nsocc = self.nelec - 2*self.ndocc
        self.maxiter = psi4.get_global_option('MAXITER')
        self.e_convergence = psi4.get_global_option('E_CONVERGENCE')
        self.nbf = mints.basisset().nbf()
        self.nsbf = 2*self.nbf
        self.vu = np.matrix(np.zeros((self.nsbf, self.nsbf)))

    # Calculate the Kronecker delta of spin-AO basis x and y

    def delta(self, x, y):

        if x < self.nbf and y < self.nbf:
            return 1
        elif not x < self.nbf and not y < self.nbf:
            return 1
        else:
            return 0                

    def build_vu(self, D):
        """
        Construct a v matrix using the density matrix and the four index integrals
        """

        r = range(int(self.nsbf))
        self.vu = np.matrix(np.zeros((self.nsbf, self.nsbf)))
        for u in r:
            u_i = u % self.nbf
            for v in r:
                v_i = v % self.nbf
                for p in r:
                    p_i = p % self.nbf
                    for q in r:
                        q_i = q % self.nbf
                        self.vu[u, v] += (self.delta(u,v) * self.delta(p,q) * self.g[u_i, p_i, v_i, q_i] - self.delta(u,q) * self.delta(p,v) * self.g[u_i, p_i, q_i, v_i]) * D[q, p] 

    def compute_energy(self):

        S = np.matrix(la.block_diag(self.S, self.S))
        T = np.matrix(la.block_diag(self.T, self.T))
        V = np.matrix(la.block_diag(self.V, self.V))
        D = np.matrix(np.zeros((self.nsbf, self.nsbf)))
        X = np.matrix(la.inv(la.sqrtm(S)))
        h = T + V
        E0 = 0
        for count in range(self.maxiter):
            F = h + self.vu
            Ft = X * F * X
            e, Ct = la.eigh(Ft)
            C = X * np.matrix(Ct)
            OC = np.matrix(C[:,:self.nelec])
            D = OC*OC.T
            self.build_vu(D)
            E1 = np.sum((np.array(h) + 0.5 * np.array(self.vu))*np.array(D.T)) + self.V_nuc
            psi4.print_out('Iteration {:<d}   {:.10f}    {:.10f}\n'.format(count, E1, E1-E0))
            if abs(E1 - E0) < self.e_convergence:
                psi4.print_out('\nFinal HF Energy: {:<5.10f}'.format(E1))
                self.C = C
                self.epsi = e
                self.ehf = E1
                break
            else:
                E0 = E1
        else:
            psi4.print_out('\n:(   Does not converge   :(')

