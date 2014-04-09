/**
 * \mainpage Basic Tensor Algebra Subroutines in C/C++ (BTAS)
 *
 * This is C++11 version of BTAS
 * Dependency to BLITZ++ library has been removed.
 *
 * T: value type
 * N: rank of array
 * Q: quantum number class
 *
 * \section FEATURES
 *
 *  - Provides generic type array classes, TArray<T, N>, STArray<T, N>, and QSTArray<T, N, Q = Quantum>
 * 
 *  - Defines DArray<N> as the alias to TArray<double, N> (using template alias in C++11). SDArray<N> and QSDArray<N, Q> as well
 * 
 *  - Provides LAPACK interfaces written in C
 * 
 *  - Provides BLAS/LAPACK-like interfaces called with DArray<N>, STArray<N>, and QSTArray<N, Q>
 * 
 *  - Provides expressive contraction, permutation, and decomposition functions
 *
 * \section COMPILATION
 *
 *  -Compiler and Library Dependencies
 * 
 * GNU GCC 4.7.0 or later
 * Intel C/C++ Compiler 13.0 or later
 *
 * BOOST library (<http://www.boost.org/>)
 * CBLAS & LAPACK library or Intel MKL library
 *
 *  - Build libbtas.a
 *
 *    cd $BTAS_ROOT/lib/
 *    make
 *
 *  - Build your code with BTAS library (GCC with MKL library)
 *
 *    g++ -std=c++0x -O3 -fopenmp -I$BTAS_ROOT/include $BTAS_ROOT/lib/libbtas.a -lboost_serialization -lmkl_core -lmkl_intel_lp64 -lmkl_sequential
 *
 * For coding, `$BTAS_ROOT/lib/tests.C` and `$BTAS_ROOT/dmrg/` involves helpful example to use BTAS
 *
 * If '-D_PRINT_WARNINGS' is specified, warning that SDArray::reserve or SDArray::insert is called with prohibited (quantum number) block is printed.
 * It gives verbose output, but helps to check undesirable behavior upon reservation and insertion.
 *
 */
#ifndef _BTAS_CXX11_QSDARRAY_H
#define _BTAS_CXX11_QSDARRAY_H 1

#include <btas/btas.h>
#include <btas/TVector.h>

#include <btas/QSPARSE/QSTArray.h>
#include <btas/QSPARSE/QSTdsum.h>

namespace btas {

//! Alias to double precision real quatum number-based sparse-array
template<size_t N, class Q = Quantum>
using QSDArray = QSTArray<double, N, Q>;

template<size_t N, class Q = Quantum>
inline void QSDdsum
(const QSDArray<N, Q>& x, const QSDArray<N, Q>& y, QSDArray<N, Q>& z) { QSTdsum(x, y, z); }

template<size_t N, size_t K, class Q = Quantum>
inline void QSDdsum
(const QSDArray<N, Q>& x, const QSDArray<N, Q>& y, const IVector<K>& trace_index, QSDArray<N, Q>& z) { QSTdsum(x, y, trace_index, z); }

}; // namespace btas

#include <btas/QSPARSE/QSDblas.h>
#include <btas/QSPARSE/QSDlapack.h>
#include <btas/QSPARSE/QSDpermute.h>
#include <btas/QSPARSE/QSDcontract.h>
//#include <btas/QSPARSE/QSDdiagonal.h>
#include <btas/QSPARSE/QSDmerge.h>

#endif // _BTAS_CXX11_QSDARRAY_H
