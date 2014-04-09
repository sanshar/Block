//
//! \file CArray.h
//! \brief Dense array for single precision complex number
//

#ifndef _BTAS_CXX11_CARRAY_H
#define _BTAS_CXX11_CARRAY_H 1

#include <iostream>
#include <iomanip>
#include <cmath>

#include <btas/TVector.h>
#include <btas/blas_cxx_interface.h>

#include <btas/DENSE/TArray.h>
#include <btas/DENSE/TSubArray.h>

namespace btas {

//! Alias to single precision complex array
template<size_t N>
using CArray = TArray<complex<float>, N>;

//! Alias to single precision complex sub-array
template<size_t N>
using CSubArray = TSubArray<complex<float>, N>;

//! Fast copying for single precision complex array (calling cblas_ccopy)
template<>
inline void _fast_copy<complex<float>>(size_t n, const complex<float>* x, complex<float>* y) {
  cblas_ccopy(n, x, 1, y, 1);
}

//! Fast scaling for single precision complex array (calling cblas_cscal)
template<>
inline void _fast_scal<complex<float>>(size_t n, const complex<float>& alpha, complex<float>* x) {
  cblas_cscal(n, &alpha, x, 1);
}

//! Fast adding  for single precision complex array (calling cblas_caxpy)
template<>
inline void _fast_add <complex<float>>(size_t n, const complex<float>* x, complex<float>* y) {
  complex<float> alpha(1.0, 0.0);
  cblas_caxpy(n, &alpha, x, 1, y, 1);
}

}; // namespace btas

//! C++ style printing function
template<size_t N>
std::ostream& operator<< (std::ostream& ost, const btas::CArray<N>& a) {
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
  for(typename btas::CArray<N>::const_iterator it = a.begin(); it != a.end(); ++it, ++n) {
    if(n % stride == 0) ost << endl << "\t";
    if(it->imag() < 0)
      ost << setw(width) << it->real() << " - " << setw(width) << fabs(it->imag()) << "i";
    else
      ost << setw(width) << it->real() << " + " << setw(width) << fabs(it->imag()) << "i";
  }
  ost << endl;

  return ost;
}

#endif // _BTAS_CXX11_CARRAY_H
