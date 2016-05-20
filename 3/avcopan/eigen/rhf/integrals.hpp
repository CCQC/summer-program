#ifndef INTEGRALS_HPP_
#define INTEGRALS_HPP_

#include <boost/shared_ptr.hpp>           // boost::shared_ptr
#include <libmints/mints.h>               // psi::OneBodyAOInt, psi::TwoBodyAOInt
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor

namespace ao_integrals {

/* typedefs */
using EigenTensor = Eigen::Tensor<double, 4>;
using EigenMatrix = Eigen::Matrix<double,-1,-1>;

/* functions */
EigenMatrix compute_oei(psi::OneBodyAOInt*);
EigenTensor compute_tei(psi::TwoBodyAOInt*);

} // end namespace ao_integrals

#endif // INTEGRALS_HPP_
