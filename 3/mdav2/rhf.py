import numpy as np
import psi4

def rhf(mol, basis='sto-3g', dE_crit=1E-6, iterlim=100):
#Takes in a psi4.core.Molecule object and a psi4 compatible basis set string.
#Returns a RHF-SCF energy and an [error/no-error] message
    
#pull in basis sets and mintshelper (molecular integrals) from psi4
    basis = psi4.core.BasisSet.build(mol, target=str(basis), key='basis')
    mints = psi4.core.MintsHelper(basis) 
    
#extract useful basic quantities from molecule, mints
    Enuc = mol.nuclear_repulsion_energy()
    natom = mol.natom()
    charge = mol.molecular_charge()
    norb = mints.basisset().nbf()
    nuclear_charges = [mol.Z(A) for A in range(natom)]
    nocc = int((sum(nuclear_charges) - charge)/2)

#Build [zero]D, overlap (S), kinetic (T), one-electron (V), two-electron (g)
#+make orthogonalizer (X) from S, make h from (T + V)
    S =  mints.ao_overlap()
    S.power(-0.5,1e-14)
    X = np.asarray(S)
    T = np.asarray(mints.ao_kinetic())
    V = np.asarray(mints.ao_potential())
    D = np.zeros(shape=T.shape)
    g = np.asarray(mints.ao_eri())
    h = T + V

#form two-e int. (g) and transpose (gt) in physicists notation
    g = g.transpose(0,2,1,3)
    gt = g.transpose(0,1,3,2)

#Begin SCF iterations. default convergence criteria is 1E-6 A.U.
    niter = 0
    E = 0
    elec_energy, dE, D, C, e = scf_iter(nocc,h,X,g,gt,D,E)
    while (niter < iterlim) and (abs(dE) > dE_crit):
        elec_energy, dE, D, C, e = scf_iter(nocc, h, X, g, gt, D, elec_energy)
        niter += 1
    E = elec_energy + Enuc
    if niter == iterlim:
        message = "Convergence Failure"
    else: 
        message = "No Convergence Failure"
    return E, message, C, g, e, nocc

def scf_iter(nocc,h,X,g,gt,D,prevE):
    nu = np.einsum('mrns,sr->mn',g - 0.5*gt,D)
    f = h + nu
    Xi = np.linalg.inv(X)
    orthog_f = X.dot(f).dot(X)
    e, eigenvectors = np.linalg.eigh(orthog_f)
    C = X.dot(eigenvectors)
    r = h.shape[0]
    
#is there a way to do this with an einsum?
    for m in range(r):
        for n in range(r):
            D[m][n] = 0
            for i in range(nocc):
                D[m][n] += C[m][i]*C[n][i]
            D[m][n] = 2*D[m][n]
    elec_energy = np.einsum('mn,nm->',(h + 0.5*nu), D)
    dE = elec_energy - prevE
    return(elec_energy, dE, D, C, e)

if __name__ == "__main__":
    mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5 
    """)
    #silence psi4 and make scf_type pk to match this code
    psi4.core.be_quiet()
    psi4.set_options({'scf_type':'pk'})

    E, msg, C, g, e, nocc = rhf(mol, basis='sto-3g', iterlim=500, dE_crit=1E-8) 
    psi_energy=psi4.energy('scf/sto-3g',molecule=mol)
    dE=E-psi_energy
    psi_match=abs(dE)<1E-6 #agreement between psi4 and this code?

    print("E: {}\nMessage: {}".format(E, msg)) 
    print("PSI4_E: {}\n".format(psi_energy))
    print("dE: {}\nMatch: {}\n".format(dE,psi_match))
