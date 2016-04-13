#ifndef RHF_HPP_
#define RHF_HPP_

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor

namespace rhf {

/* typedefs */
template<class Class>
using Shared   = boost::shared_ptr<Class>;
using Molecule = psi::Molecule;
using Options  = psi::Options;
using Mints    = psi::MintsHelper;

using EigenTensor  = Eigen::Tensor<double, 4>;
using EigenMatrix  = Eigen::Matrix<double,-1,-1>;
using Eigenvalues  = Eigen::Matrix<double,-1, 1>;

class RHF {
 private:
  /* class variables */
  int norb_;
  int nocc_;
  double Vnu_;
  Options options_;
  EigenMatrix h_; // core hamiltonian
  EigenMatrix X_; // orthogonalizer
  EigenTensor g_; // two-electron integrals in physicist's notation, <mu nu | rh si>

 public:
  /* constructor */
  RHF(Shared<Molecule> mol, Options options);
  /* class methods */
  double compute_energy();
};

} // end namespace rhf

#endif
