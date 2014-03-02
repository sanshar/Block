//
//! \file DArray.h
//! \brief Dense array for double precision real number
//

#ifndef _BTAS_CXX11_DARRAY_H
#define _BTAS_CXX11_DARRAY_H 1

#include <iostream>
#include <iomanip>

#include <btas/TVector.h>
#include <btas/blas_cxx_interface.h>

#include <btas/DENSE/TArray.h>
#include <btas/DENSE/TSubArray.h>

namespace btas {

//! Alias to double precision real array
template<size_t N>
using DArray = TArray<double, N>;

//! Alias to double precision real sub-array
template<size_t N>
using DSubArray = TSubArray<double, N>;

//! Fast copying for double precision real array (calling cblas_dcopy)
template<>
inline void _fast_copy<double>(size_t n, const double* x, double* y) { cblas_dcopy(n, x, 1, y, 1); }

//! Fast scaling for double precision real array (calling cblas_dscal)
template<>
inline void _fast_scal<double>(size_t n, const double& alpha, double* x) { cblas_dscal(n, alpha, x, 1); }

//! Fast adding  for double precision real array (calling cblas_daxpy)
template<>
inline void _fast_add <double>(size_t n, const double* x, double* y) { cblas_daxpy(n, 1.0, x, 1, y, 1); }

}; // namespace btas

//! C++ style printing function
template<size_t N>
std::ostream& operator<< (std::ostream& ost, const btas::DArray<N>& a) {
  using std::setw;
  using std::endl;
  // detect ostream status for floating point value
  int width = ost.precision() + 4;
  if(ost.flags() & std::ios::scientific)
    width += 4;
  else
    ost.setf(std::ios::fixed, std::ios::floatfield);

  // printing array shape
  const btas::IVector<N>& a_shape = a.shape();
  ost << "shape [ ";
  for(int i = 0; i < N-1; ++i) ost << a_shape[i] << " x ";
  ost << a_shape[N-1] << " ] " << endl;
  ost << "----------------------------------------------------------------------------------------------------" << endl;
  // printing array elements
  int stride = a.shape(N-1);
  int n = 0;
  for(typename btas::DArray<N>::const_iterator it = a.begin(); it != a.end(); ++it, ++n) {
    if(n % stride == 0) ost << endl << "\t";
    ost << setw(width) << *it;
  }
  ost << endl;

  return ost;
}

#include <btas/DENSE/Dblas.h>
#include <btas/DENSE/Dlapack.h>
#include <btas/DENSE/Dpermute.h>
#include <btas/DENSE/Dcontract.h>
#include <btas/DENSE/Ddiagonal.h>

#endif // _BTAS_CXX11_DARRAY_H
