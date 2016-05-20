#include "rhf.hpp"
#include "integrals.hpp"
#include <stdexcept> // std::invalid_argument


namespace rhf {

/* constructor */
RHF::RHF(Shared<Molecule> mol, Options options) : options_(options) {
  Shared<Mints> mints(new Mints(options));

  Vnu_  = mol->nuclear_repulsion_energy();
  norb_ = mints->basisset()->nbf();

  int nelec = - mol->molecular_charge();
  for (int A = 0; A < mol->natom(); ++A) nelec += mol->Z(A);
  nocc_ = nelec / 2;

  EigenMatrix S = ao_integrals::compute_oei(mints->integral()->ao_overlap()  );
  EigenMatrix T = ao_integrals::compute_oei(mints->integral()->ao_kinetic()  );
  EigenMatrix V = ao_integrals::compute_oei(mints->integral()->ao_potential());
  EigenTensor G = ao_integrals::compute_tei(mints->integral()->eri()         );

  // build core Hamiltonian
  h_ = T + V;

  // build orthogonalizer
  Eigen::SelfAdjointEigenSolver<EigenMatrix> eigS(S);
  Eigenvalues x = eigS.eigenvalues().array().pow(-1./2);
  EigenMatrix U = eigS.eigenvectors();
  X_ = U * x.asDiagonal() * U.transpose();

  // get physicist's notation integrals
  std::array<int,4> phys({0,2,1,3});
  g_ = G.shuffle(phys);

  // if it isn't a closed-shell molecule, complain
  if (nelec % 2 != 0) throw std::invalid_argument("RHF code requires a closed-shell molecule.");
}


/* methods */
double RHF::compute_energy() {
  std::cout << "Hello, Brave New World!" << std::endl;

  EigenMatrix D = EigenMatrix::Zero(norb_, norb_);

  for(int iter = 0; iter < options_.get_int("MAXITER"); ++iter) {
    std::cout << iter << std::endl;
  }

  return 0.0;
}


} // end namespace rhf
