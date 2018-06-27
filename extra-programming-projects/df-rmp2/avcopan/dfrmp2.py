import psi4
import numpy
import scipy.linalg

BASIS = 'aug-cc-pvdz'
AUXBASIS = 'aug-cc-pvdz-ri'
CHARGE = 0
SPIN = 0
LABELS = ('O', 'H', 'H')
COORDS = ((0.000000000000,  0.000000000000, -0.143225816552),
          (0.000000000000,  1.638036840407,  1.136548822547),
          (0.000000000000, -1.638036840407,  1.136548822547))


def main():
    c, e, na = restricted_orbitals(basis=BASIS, labels=LABELS, coords=COORDS,
                                   charge=CHARGE, spin=SPIN)
    fr = factorized_repulsion_integrals(basis=BASIS, auxbasis=AUXBASIS,
                                        labels=LABELS, coords=COORDS)

    x = numpy.newaxis
    eo, ev = numpy.split(e, (na,))
    co, cv = numpy.split(c, (na,), axis=1)

    riax = numpy.tensordot(fr, co, axes=(0, 0))
    riax = numpy.tensordot(riax, cv, axes=(0, 0))
    riax = numpy.transpose(riax, (1, 2, 0))

    riajb = numpy.tensordot(riax, riax, axes=(-1, -1))
    gijab = numpy.transpose(riajb, (0, 2, 1, 3))
    gijba = numpy.transpose(riajb, (0, 2, 3, 1))
    t2 = (2*gijab - gijba) / (+ eo[:, x, x, x] + eo[x, :, x, x]
                              - ev[x, x, :, x] - ev[x, x, x, :])

    e_corr = numpy.vdot(t2, gijab)
    print(e_corr)


def restricted_orbitals(basis, labels, coords, charge=0, spin=0, niter=100,
                        e_thresh=1e-12, d_thresh=1e-9, guess='gwh'):
    """restricted Hartree-Fock orbitals

    :param basis: basis set name
    :type basis: str
    :param labels: atomic symbols labeling the nuclei
    :type labels: tuple
    :param coords: nuclear coordinates in Bohr
    :type coords: numpy.ndarray
    :param charge: total molecular charge
    :type charge: int
    :param spin: number of unpaired electrons
    :type spin: int
    :param niter: maximum number of iterations
    :type niter: int
    :param e_thresh: energy convergence threshold
    :type e_thresh: float
    :param d_thresh: density convergence threshold
    :type d_thresh: float
    :param guess: hartree-fock starting guess
    :type guess: str

    :return: orbital coefficients and energies, along with the occupation count
    :rtype: (numpy.ndarray, numpy.ndarray, int)
    """
    psi4.set_options({'e_convergence': e_thresh, 'd_convergence': d_thresh,
                      'maxiter': niter, 'guess': guess, 'reference': 'RHF'})
    mol = psi4_molecule(labels=labels, coords=coords, charge=charge, spin=spin)
    wfn = psi4.core.Wavefunction.build(mol, basis)
    sf, _ = psi4.driver.dft_funcs.build_superfunctional("HF", False)
    hf = psi4.core.RHF(wfn, sf)
    hf.compute_energy()
    return numpy.array(hf.Ca()), numpy.array(hf.epsilon_a()), int(hf.nalpha())


def factorized_repulsion_integrals(basis, auxbasis, labels, coords):
    """RI-factorized electron-electron repulsion integrals

    :param basis: basis set name
    :type basis: str
    :param basis: auxiliary basis set name
    :type basis: str
    :param labels: atomic symbols labeling the nuclei
    :type labels: tuple
    :param coords: nuclear coordinates in Bohr
    :type coords: numpy.ndarray

    :return: a three-index tensor; the last axis is the auxiliary one
    :rtype: numpy.ndarray
    """
    j = coulomb_metric_integrals(basis=auxbasis, labels=labels, coords=coords)
    r3c = threecenter_repulsion_integrals(
            basis=basis, auxbasis=auxbasis, labels=labels, coords=coords)

    x = scipy.linalg.cholesky(scipy.linalg.inv(j))
    return numpy.tensordot(r3c, x, axes=(2, 1))


def threecenter_repulsion_integrals(basis, auxbasis, labels, coords):
    """three-center electron-electron repulsion integrals

    :param basis: basis set name
    :type basis: str
    :param basis: auxiliary basis set name
    :type basis: str
    :param labels: atomic symbols labeling the nuclei
    :type labels: tuple
    :param coords: nuclear coordinates in Bohr
    :type coords: numpy.ndarray

    :return: a three-index tensor; the last axis is the auxiliary one
    :rtype: numpy.ndarray
    """
    psi4.core.clean()
    mol = psi4_molecule(labels=labels, coords=coords)
    bs = psi4.core.BasisSet.build(mol, 'BASIS', basis)
    bs_aux = psi4.core.BasisSet.build(mol, 'DF_BASIS_MP2', auxbasis)
    bs_empty = psi4.core.BasisSet.zero_ao_basis_set()
    mh = psi4.core.MintsHelper(bs)
    return numpy.squeeze(mh.ao_eri(bs, bs, bs_aux, bs_empty))


def coulomb_metric_integrals(basis, labels, coords):
    """coulomb metric integrals

    :param basis: basis set name
    :type basis: str
    :param labels: atomic symbols labeling the nuclei
    :type labels: tuple
    :param coords: nuclear coordinates in Bohr
    :type coords: numpy.ndarray

    :return: a four-index tensor of equal dimensions
    :rtype: numpy.ndarray
    """
    psi4.core.clean()
    mol = psi4_molecule(labels=labels, coords=coords)
    bs = psi4.core.BasisSet.build(mol, 'BASIS', basis)
    bs_empty = psi4.core.BasisSet.zero_ao_basis_set()
    mh = psi4.core.MintsHelper(bs)
    return numpy.squeeze(mh.ao_eri(bs, bs_empty, bs, bs_empty))


def psi4_molecule(labels, coords, charge=0, spin=0):
    """
    :param labels: atomic symbols labeling the nuclei
    :type labels: tuple
    :param coords: nuclear coordinates in Bohr
    :type coords: numpy.ndarray
    :param charge: total molecular charge
    :type charge: int
    :param spin: number of unpaired electrons
    :type spin: int

    :rtype: psi4.core.Molecule
    """
    xyz = coordinate_string(labels=labels, coords=coords)
    mol = psi4.core.Molecule.create_molecule_from_string(xyz)
    mol.set_molecular_charge(charge)
    mol.set_multiplicity(spin + 1)
    mol.reset_point_group('C1')
    mol.update_geometry()
    return mol


def coordinate_string(labels, coords):
    """coordinate string

    :param labels: atomic symbols labeling the nuclei
    :type labels: tuple
    :param coords: nuclear coordinates in Bohr
    :type coords: numpy.ndarray

    :returns: the molecular geometry
    :rtype: str
    """
    line = '{:2s} {: >17.12f} {: >17.12f} {: >17.12f}'
    xyz = '\n'.join(line.format(l, *c) for l, c in zip(labels, coords))
    return 'units bohr\n{:s}'.format(xyz)


if __name__ == '__main__':
    main()
