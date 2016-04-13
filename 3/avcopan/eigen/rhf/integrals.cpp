#include "integrals.hpp"

namespace ao_integrals {

EigenMatrix compute_oei(psi::OneBodyAOInt* O) {
  boost::shared_ptr<psi::BasisSet> bs1 = O->basis1(); int dim1 = bs1->nbf();
  boost::shared_ptr<psi::BasisSet> bs2 = O->basis2(); int dim2 = bs2->nbf();
 
  EigenMatrix A(dim1, dim2);
  double* pA = A.data();

  const double* buffer = O->buffer();

  for (int M = 0; M < bs1->nshell(); ++M) {
    for (int N = 0; N < bs2->nshell(); ++N) {

      O->compute_shell(M, N);

      int index = 0;
      for (int m = 0; m < bs1->shell(M).nfunction(); ++m)
        for (int n = 0; n < bs2->shell(N).nfunction(); ++n)
          pA[(bs1->shell(M).function_index() + m)*1   +
             (bs2->shell(N).function_index() + n)*dim1] = buffer[index++];
    }
  }

  return A;
}

EigenTensor compute_tei(psi::TwoBodyAOInt* T) {

  boost::shared_ptr<psi::BasisSet> bs1 = T->basis1(); int dim1 = bs1->nbf();
  boost::shared_ptr<psi::BasisSet> bs2 = T->basis2(); int dim2 = bs2->nbf();
  boost::shared_ptr<psi::BasisSet> bs3 = T->basis3(); int dim3 = bs3->nbf();
  boost::shared_ptr<psi::BasisSet> bs4 = T->basis4(); int dim4 = bs4->nbf();

  EigenTensor A(dim1, dim2, dim3, dim4);
  double* pA = A.data();

  const double* buffer = T->buffer();

  for (int M = 0; M < bs1->nshell(); ++M) {
    for (int N = 0; N < bs2->nshell(); ++N) {
      for (int P = 0; P < bs3->nshell(); ++P) {
        for (int Q = 0; Q < bs4->nshell(); ++Q) {

          T->compute_shell(M,N,P,Q);

          int index = 0;
          for (int m = 0; m < bs1->shell(M).nfunction(); ++m)
            for (int n = 0; n < bs2->shell(N).nfunction(); ++n)
              for (int p = 0; p < bs3->shell(P).nfunction(); ++p)
                for (int q = 0; q < bs4->shell(Q).nfunction(); ++q)
                  pA[(bs1->shell(M).function_index() + m)*1             +
                     (bs2->shell(N).function_index() + n)*dim1          +
                     (bs3->shell(P).function_index() + p)*dim1*dim2     +
                     (bs4->shell(Q).function_index() + q)*dim1*dim2*dim3] = buffer[index++];

        }
      }
    }
  }

  return A;
}

} // end namespace ao_integrals
