# Alternate UHF code separating alpha and beta spin orbitals

import psi4
import numpy as np
import configparser 
import scipy.linalg as la

class UHF:
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
        self.nelec = -self.mol.molecular_charge()
        for A in range(self.mol.natom()):
            self.nelec += int(self.mol.Z(A))
        self.nocc = self.nelec 
        self.ntot = self.mints.basisset().nbf()
        self.S = self.mints.ao_overlap().to_array()
        self.T = self.mints.ao_kinetic().to_array()
        self.V = self.mints.ao_potential().to_array()
        self.I = self.mints.ao_eri().to_array() 
        
        self.H = self.T + self.V
        self.C = np.zeros_like((len(self.H)*2,len(self.H)*2))
        self.C_a = np.zeros_like(self.H)
        self.C_b = np.zeros_like(self.H)
        self.ea = np.zeros(len(self.H))
        self.eb = np.zeros(len(self.H))
        self.e = np.zeros(len(self.H)*2)
        A = self.mints.ao_overlap()
        A.power(-0.5, 1.e-16)
        self.A = A.to_array()
        
    def get_energy(self):
        mol, max_iter, nalpha,nbeta, S, T, V, I, H, A, C_a, C_b, ea, eb = self.mol, self.max_iter, self.nalpha, self.nbeta, self.S, self.T, self.V, self.I, self.H, self.A, self.C_a, self.C_b, self.ea, self.eb
        
        Fa = H
        Fb = H 
        E_old = 0.0
        Da_old = np.zeros_like(H)
        Db_old = np.zeros_like(H)
        for iteration in range(1, self.max_iter+1):


            Ft_a = A.dot(Fa).dot(A)
            Ft_b = A.dot(Fb).dot(A)
            ea, C_a = np.linalg.eigh(Ft_a)
            eb, C_b = np.linalg.eigh(Ft_b)
            C_a = A.dot(C_a)
            C_b = A.dot(C_b)
            Ca = C_a[:,:nalpha]
            Cb = C_b[:,:nbeta]
            Da = np.einsum('pi, qi->pq', Ca, Ca)
            Db = np.einsum('pi, qi->pq', Cb, Cb)
            
            Ja = np.einsum('pqrs, rs->pq', I, Da)
            Jb = np.einsum('pqrs, rs->pq', I, Db)
            Ka = np.einsum('prqs, rs->pq', I, Da)
            Kb = np.einsum('prqs, rs->pq', I, Db)
            Fa = H + Ja - Ka + Jb
            Fb = H + Jb - Kb + Ja
        
            E_SCF = (1/2)*(np.einsum('pq, pq->', Fa+H, Da) + np.einsum('pq,pq->',Fb+H, Db)) +mol.nuclear_repulsion_energy()
            print('UHF iteration {:3d}: energy {:20.14f} dE {:1.5E}'.format(iteration, E_SCF, (E_SCF - E_old)))

            if (abs(E_SCF - E_old) < 1.e-10):
                break
             
            E_old = E_SCF
            Da_old = Da
            Db_old = Db

        self.C_a, self.C_b = C_a, C_b 
        self.ea, self.eb = ea, eb
