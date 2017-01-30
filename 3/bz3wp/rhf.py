import psi4
import numpy as np
import configparser 
import scipy.linalg as la

class RHF:
    def __init__(self, filename = 'Options.ini'):
        config = configparser.ConfigParser()
        config.read(filename)
        self.mol = psi4.geometry(config['DEFAULT']['molecule'])
        self.mol.update_geometry()
        basis = psi4.core.BasisSet.build(self.mol, 'BASIS', config['DEFAULT']['basis'])
        mints = psi4.core.MintsHelper(basis)
        self.max_iter = int(config['SCF']['max_iter'])
        self.ntot = mints.basisset().nbf() 
        self.nelec = -self.mol.molecular_charge()
        for A in range(self.mol.natom()):
            self.nelec += int(self.mol.Z(A))
        self.ndocc = int(self.nelec/2)

        S = mints.ao_overlap().to_array()
        T = mints.ao_kinetic().to_array()
        V = mints.ao_potential().to_array()
        self.G = mints.ao_eri().to_array() 
        
        self.H = T + V
        self.C = np.zeros_like(self.H)
        self.e = np.zeros(len(self.H))
        self.E_SCF = 0.0 
        self.A = np.matrix(la.inv(la.sqrtm(S)))
        
    def get_energy(self):
        mol, max_iter, ndocc, H, G, A, C, e = self.mol, self.max_iter, self.ndocc,self.H, self.G, self.A, self.C, self.e
        E_old = 0.0
        D_old = np.zeros_like(H)

        for iteration in range(1, self.max_iter+1):
            
            J = np.einsum('pqrs, rs->pq', G, D_old)
            K = np.einsum('prqs, rs->pq', G, D_old)
            F = H + J*2 - K
            
            Ft = A.dot(F).dot(A)
            e, C = np.linalg.eigh(Ft)
            C = A.dot(C)
            Cocc = C[:,:ndocc]
            D = np.einsum('pi, qi->pq', Cocc, Cocc)
            
            E_SCF = np.einsum('pq, pq->', F+H, D) + mol.nuclear_repulsion_energy()
            D_norm = np.linalg.norm(D-D_old)
            print('RHF iteration {:3d}: energy {:20.14f} dE {:2.5E} D_norm {:2.5E}'.format(iteration, E_SCF, (E_SCF - E_old), D_norm))

            if (abs(E_SCF - E_old) < 1.e-10) and D_norm < 1.e-10:
                break
             
            E_old = E_SCF
            D_old = D
        
        self.e = e
        self.C = C
        self.E_SCF = E_SCF    
