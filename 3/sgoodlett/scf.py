from geom import molecule
import numpy as np
from masses import get_mass, get_charge
import psi4
from scipy.linalg import sqrtm, inv, eigh

class SCF():

    def __init__(self,mol,basis,convergence=6,maxiter=50):
        self.basis = basis
        self.convergence = convergence
        self.maxiter = maxiter
        psi4.set_memory('500 MB')
        psi4.set_output_file('output_scf.dat',False)
        self.mol = psi4.geometry(mol.xyz_string(option=1))
        psi4.set_options({'basis':self.basis,'scf_type':'pk'})
        #psi4.energy('scf')
        self.wfn = psi4.core.Wavefunction.build(self.mol,psi4.core.get_global_option('BASIS'))
        self.mints = psi4.core.MintsHelper(self.wfn.basisset())

        #Init integrals
        self.S = np.matrix(self.mints.ao_overlap())
        self.T = np.matrix(self.mints.ao_kinetic())
        self.V = np.matrix(self.mints.ao_potential())
        self.g_chem = np.array(self.mints.ao_eri())
        self.g = self.g_chem.transpose((0,2,1,3))
        self.gt = self.g.transpose((0,1,3,2))
        self.X = np.matrix(inv(sqrtm(self.S)))

        #Init Nuclear Repulsion Energy
        self.E_nuc = self.mol.nuclear_repulsion_energy()
        self.natom = self.mol.natom()
        self.charge = self.mol.molecular_charge()
        self.norbitals = self.mints.basisset().nbf()

        self.n_e = sum(mol.charges) - self.charge
        try:
            self.n_occ = int(self.n_e / 2)
        except:
            print('Error: Not all orbitals fully occupied, try using UHF instead')
            pass
        self.D = np.zeros((self.norbitals,self.norbitals)) #Might be wrong???

    def scf_cycle(self):

        #Build Fock matrix (eq 4)
        v = np.einsum('ijkl,lj->ik',self.g - 0.5 * self.gt, self.D)
        h = self.T + self.V
        f = h + v

        #Transform f -> ~f to orthogonalized AO basis (eq 6)
        f_ortho = np.matmul(self.X,np.matmul(f,self.X))

        #Diagonalize ~f, yielding orbital energies \epsilon_p and MO coeffecients ~C_{\mu p} (eq. 6)
        eps, C_ortho = eigh(f_ortho)

        #Backtransform ~C -> C to original AO basis (eq. 6)
        C = np.matmul(self.X,np.matrix(C_ortho))

        #Build Density matrix D (eq. 3)
        self.D = 2 * np.einsum('ij,kj->ik',C[:,:self.n_occ],np.conj(C[:,:self.n_occ]))

        #Evaluate Energy (eq. 5)
        hv = h + (0.5 * v)
        energy = np.einsum('ij,ji->',hv,self.D) + self.E_nuc

        return energy, C, eps

    def energy(self,return_integrals=False):
        index = 0
        old_energy = 0
        eps = 10 ** (-1 * self.convergence)
        while True:
            energy, C, orbital_e = self.scf_cycle()

            if abs(energy - old_energy) < eps:
                if return_integrals:
                    return energy, C, orbital_e, self.g_chem
                else:
                    return energy, C, orbital_e
            if index > self.maxiter:
                print('Surpassed SCF max iter')
                break

            old_energy = energy
            index += 1


if __name__ == '__main__':
    mol = molecule(units='Bohr')
    mol.to_angstrom()
    e = SCF(mol,'cc-pVDZ',convergence=10,maxiter=1000)
    e.energy()