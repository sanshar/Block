//
//! \file SArray.h
//! \brief Dense array for single precision real number
//

#ifndef _BTAS_CXX11_SARRAY_H
#define _BTAS_CXX11_SARRAY_H 1

#include <iostream>
#include <iomanip>

#include <btas/TVector.h>
#include <btas/blas_cxx_interface.h>

#include <btas/DENSE/TArray.h>
#include <btas/DENSE/TSubArray.h>

namespace btas {

//! Alias to single precision real array
template<size_t N>
using SArray = TArray<float, N>;

//! Alias to single precision real sub-array
template<size_t N>
using SSubArray = TSubArray<float, N>;

//! Fast copying for single precision real array (calling cblas_scopy)
template<>
inline void _fast_copy<float>(size_t n, const float* x, float* y) { cblas_scopy(n, x, 1, y, 1); }

//! Fast scaling for single precision real array (calling cblas_sscal)
template<>
inline void _fast_scal<float>(size_t n, const float& alpha, float* x) { cblas_sscal(n, alpha, x, 1); }

//! Fast adding  for single precision real array (calling cblas_saxpy)
template<>
inline void _fast_add <float>(size_t n, const float* x, float* y) { cblas_saxpy(n, 1.0, x, 1, y, 1); }

}; // namespace btas

//! C++ style printing function
template<size_t N>
std::ostream& operator<< (std::ostream& ost, const btas::SArray<N>& a) {
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
  for(typename btas::SArray<N>::const_iterator it = a.begin(); it != a.end(); ++it, ++n) {
    if(n % stride == 0) ost << endl << "\t";
    ost << setw(width) << *it;
  }
  ost << endl;

  return ost;
}

#endif // _BTAS_CXX11_SARRAY_H
