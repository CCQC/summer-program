import psi4
import numpy as np
import configparser 

class RHF:
    def __init__(self, filename = 'Options.ini'):
        self.config = configparser.ConfigParser()
        self.config.read(filename)
        self.mol = psi4.geometry(self.config['DEFAULT']['molecule'])
        self.mol.update_geometry()
        self.basis = psi4.core.BasisSet.build(self.mol, 'BASIS', self.config['DEFAULT']['basis'])
        self.mints = psi4.core.MintsHelper(self.basis)
        self.max_iter = int(self.config['SCF']['max_iter'])
        self.nalpha = int(self.config['DEFAULT']['nalpha'])
        self.nbeta = int(self.config['DEFAULT']['nbeta'])
        self.ntot = self.mints.basisset().nbf() 
        if (self.nalpha + self.nbeta)/2 % 1:
            print('RHF not sufficient for open shell system')
        else:
            self.ndocc = int((self.nalpha + self.nbeta)/2)

        self.S = self.mints.ao_overlap().to_array()
        self.T = self.mints.ao_kinetic().to_array()
        self.V = self.mints.ao_potential().to_array()
        self.I = self.mints.ao_eri().to_array() 
        
        self.H = self.T + self.V
        self.C = np.zeros_like(self.H)
        self.e = np.zeros(len(self.H))
        A = self.mints.ao_overlap()
        A.power(-0.5, 1.e-16)
        self.A = A.to_array()
        
    def get_energy(self):
        mol, max_iter, ndocc, S, T, V, I, H, A, C, e = self.mol, self.max_iter, self.ndocc, self.S, self.T, self.V, self.I, self.H, self.A, self.C, self.e
        
        F = H 
        E_old = 0.0
        D_old = np.zeros_like(H)
        for iteration in range(1, self.max_iter+1):


            Ft = A.dot(F).dot(A)
            e, C = np.linalg.eigh(Ft)
            C = A.dot(C)
            Cocc = C[:,:ndocc]
            D = np.einsum('pi, qi->pq', Cocc, Cocc)
            
            J = np.einsum('pqrs, rs->pq', I, D)
            K = np.einsum('prqs, rs->pq', I, D)
            F = H + J*2 - K
        
            E_SCF = np.einsum('pq, pq->', F+H, D) + mol.nuclear_repulsion_energy()
            D_norm = np.linalg.norm(D-D_old)
            print('RHF iteration {:3d}: energy {:20.14f} dE {:2.5E} D_norm {:2.5E}'.format(iteration, E_SCF, (E_SCF - E_old), D_norm))

            if (abs(E_SCF - E_old) < 1.e-10) and D_norm < 1.e-10:
                break
             
            E_old = E_SCF
            D_old = D
        
        self.e = e
        self.C = C
            
