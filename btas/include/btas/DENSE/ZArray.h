//
//! \file ZArray.h
//! \brief Dense array for double precision complex number
//

#ifndef _BTAS_CXX11_ZARRAY_H
#define _BTAS_CXX11_ZARRAY_H 1

#include <iostream>
#include <iomanip>
#include <cmath>

#include <btas/TVector.h>
#include <btas/blas_cxx_interface.h>

#include <btas/DENSE/TArray.h>
#include <btas/DENSE/TSubArray.h>

namespace btas {

//! Alias to double precision complex array
template<size_t N>
using ZArray = TArray<complex<double>, N>;

//! Alias to double precision complex sub-array
template<size_t N>
using ZSubArray = TSubArray<complex<double>, N>;

//! Fast copying for double precision complex array (calling cblas_zcopy)
template<>
inline void _fast_copy<complex<double>>(size_t n, const complex<double>* x, complex<double>* y) {
 cblas_zcopy(n, x, 1, y, 1);
}

//! Fast scaling for double precision complex array (calling cblas_zscal)
template<>
inline void _fast_scal<complex<double>>(size_t n, const complex<double>& alpha, complex<double>* x) {
  cblas_zscal(n, &alpha, x, 1);
}

//! Fast adding  for double precision complex array (calling cblas_zaxpy)
template<>
inline void _fast_add <complex<double>>(size_t n, const complex<double>* x, complex<double>* y) {
  complex<double> alpha(1.0, 0.0);
  cblas_zaxpy(n, &alpha, x, 1, y, 1);
}

}; // namespace btas

//! C++ style printing function
template<size_t N>
std::ostream& operator<< (std::ostream& ost, const btas::ZArray<N>& a) {
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
  for(typename btas::ZArray<N>::const_iterator it = a.begin(); it != a.end(); ++it, ++n) {
    if(n % stride == 0) ost << endl << "\t";
    if(it->imag() < 0)
      ost << setw(width) << it->real() << " - " << setw(width) << fabs(it->imag()) << "i";
    else
      ost << setw(width) << it->real() << " + " << setw(width) << fabs(it->imag()) << "i";
  }
  ost << endl;

  return ost;
}

#endif // _BTAS_CXX11_ZARRAY_H
