import psi4
import numpy as np
import configparser 
import scipy.linalg as la

class UHF:
    def __init__(self, filename = 'Options.ini'):
        config = configparser.ConfigParser()
        config.read(filename)
        self.mol = psi4.geometry(config['DEFAULT']['molecule'])
        self.mol.update_geometry()
        basis = psi4.core.BasisSet.build(self.mol,'BASIS', config['DEFAULT']['basis'])
        mints = psi4.core.MintsHelper(basis)
        self.max_iter = int(config['SCF']['max_iter'])
        self.nelec = -self.mol.molecular_charge()
        for A in range(self.mol.natom()):
            self.nelec += int(self.mol.Z(A))
        self.nocc = self.nelec
        self.nto = 2* mints.basisset().nbf()
        
        S = block_oei(mints.ao_overlap())
        T = block_oei(mints.ao_kinetic().to_array())
        V = block_oei(mints.ao_potential().to_array())
        G = block_tei(mints.ao_eri().to_array())
        
        
        self.H = T + V
        self.g = G.transpose(0,2,1,3) - G.transpose(0,2,3,1)
        
        self.C = np.zeros_like(self.H)
        self.e = np.zeros(len(self.H))

        self. A = np.matrix(la.inv(la.sqrtm(S)))

    def get_energy(self):
        mol, max_iter, g, H, A, C, e, nocc = self.mol, self.max_iter,self.g, self.H, self.A, self.C, self.e, self.nocc
        
        E_old = 0.0
        D_old = np.zeros_like(H)
        
        for iteration in range(1, self.max_iter+1):
            
            Gao = np.einsum('pqrs,sq', g, D_old)
            F = H + Gao
            Ft = A.dot(F).dot(A)
            e, C = np.linalg.eigh(Ft)
            C = A.dot(C)
            Cnocc = C[:,:nocc]
            D = np.einsum('pi, qi->pq', Cnocc, Cnocc)
            
            E_SCF = (1/2)*(np.einsum('pq, pq->', F+H, D)) + mol.nuclear_repulsion_energy()
            print('UHF iteration {:3d}: energy {:20.14f} dE {:1.5E}'.format(iteration, E_SCF, (E_SCF - E_old)))

            if (abs(E_SCF - E_old) < 1.e-10):
                break
             
            E_old = E_SCF
            D_old = D

        self.e = e
        self.C = C

def block_oei(A):
    A = la.block_diag(A, A)
    return A

def block_tei(A):
    I = np.identity(2)
    A = np.kron(I, A)
    return np.kron(I, A.T)


