//
/*! \file TVector.h
 *  \brief Utility functions for fixed and variable size vectors.
 *
 *  Boost array or STL array (C++11) class is used as a fixed size vector TVector
 *  Convenient shaping function, dot product, sequence constructor, transposition,
 *  and permutation functions are provided.
 *
 *  C/C++ interfaces of fast copying, scaling, and adding are provided, which are
 *  specialized for each variable type to call corresponding BLAS function.
 */

#ifndef _BTAS_CXX11_TVECTOR_H
#define _BTAS_CXX11_TVECTOR_H 1

#include <vector>
#include <array>
#include <complex>
#include <algorithm>

#include <boost/serialization/serialization.hpp>

#include <btas/btas.h>

namespace boost {
namespace serialization {

//####################################################################################################
// Boost serialization for std::array
//####################################################################################################

//! Enables to use boost serialization
template<class Archive, typename T, size_t N>
void serialize(Archive& ar, std::array<T, N>& vec, const unsigned int version) {
  for(size_t i = 0; i < N; ++i) ar & vec[i];
}

}; // namespace serialization
}; // namespace boost

namespace btas {

//####################################################################################################
// Template aliases to fixed-rank array
//####################################################################################################

//! Template aliases to std::array<T, N>, for convenience
template<typename T, size_t N>
using TVector = std::array<T, N>;

//! Template aliases to std::array<int, N>, for convenience
template<size_t N>
using IVector = std::array<int, N>;

//! Convenient constructor of TVector with const value
template<typename T, size_t N>
TVector<T, N> uniform(const T& value) {
  TVector<T, N> vec; vec.fill(value);
  return std::move(vec);
}

//####################################################################################################
// Convenient TVector constructor: make_array
//####################################################################################################

//! Convenient TVector constructor for N = 1
template<typename T>
inline TVector<T, 1> make_array(T v01) {
  TVector<T, 1> _tvec = { v01 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 2
template<typename T>
inline TVector<T, 2> make_array(T v01, T v02) {
  TVector<T, 2> _tvec = { v01, v02 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 3
template<typename T>
inline TVector<T, 3> make_array(T v01, T v02, T v03) {
  TVector<T, 3> _tvec = { v01, v02, v03 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 4
template<typename T>
inline TVector<T, 4> make_array(T v01, T v02, T v03, T v04) {
  TVector<T, 4> _tvec = { v01, v02, v03, v04 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 5
template<typename T>
inline TVector<T, 5> make_array(T v01, T v02, T v03, T v04, T v05) {
  TVector<T, 5> _tvec = { v01, v02, v03, v04, v05 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 6
template<typename T>
inline TVector<T, 6> make_array(T v01, T v02, T v03, T v04, T v05, T v06) {
  TVector<T, 6> _tvec = { v01, v02, v03, v04, v05, v06 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 7
template<typename T>
inline TVector<T, 7> make_array(T v01, T v02, T v03, T v04, T v05, T v06, T v07) {
  TVector<T, 7> _tvec = { v01, v02, v03, v04, v05, v06, v07 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 8
template<typename T>
inline TVector<T, 8> make_array(T v01, T v02, T v03, T v04, T v05, T v06, T v07, T v08) {
  TVector<T, 8> _tvec = { v01, v02, v03, v04, v05, v06, v07, v08 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 9
template<typename T>
inline TVector<T, 9> make_array(T v01, T v02, T v03, T v04, T v05, T v06, T v07, T v08, T v09) {
  TVector<T, 9> _tvec = { v01, v02, v03, v04, v05, v06, v07, v08, v09 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 10
template<typename T>
inline TVector<T,10> make_array(T v01, T v02, T v03, T v04, T v05, T v06, T v07, T v08, T v09, T v10) {
  TVector<T,10> _tvec = { v01, v02, v03, v04, v05, v06, v07, v08, v09, v10 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 11
template<typename T>
inline TVector<T,11> make_array(T v01, T v02, T v03, T v04, T v05, T v06, T v07, T v08, T v09, T v10, T v11) {
  TVector<T,11> _tvec = { v01, v02, v03, v04, v05, v06, v07, v08, v09, v10, v11 };
  return std::move(_tvec);
}

//! Convenient TVector constructor for N = 12
template<typename T>
inline TVector<T,12> make_array(T v01, T v02, T v03, T v04, T v05, T v06, T v07, T v08, T v09, T v10, T v11, T v12) {
  TVector<T,12> _tvec = { v01, v02, v03, v04, v05, v06, v07, v08, v09, v10, v11, v12 };
  return std::move(_tvec);
}

//####################################################################################################
// Overloaded operator for index product of TVector<std::vector<T>, N>
//####################################################################################################

//! Indexed product of TVector<std::vector<T>, N>
/*! e.g. index = { i, j, k, l }
 *  return vec[0][i] * vec[1][j] * vec[2][k] * vec[3][l] */
template<typename T, size_t N>
inline T operator* (const TVector<std::vector<T>, N>& vec, const IVector<N>& index) {
  T prod = vec[0][index[0]];
  for(int i = 1; i < N; ++i) prod = prod * vec[i][index[i]];
  return prod;
}

//! Indexed product of TVector<std::vector<T>, N>
/*! e.g. index = { i, j, k, l }
 *  return vec[0][i] * vec[1][j] * vec[2][k] * vec[3][l] */
template<typename T, size_t N>
inline T operator* (const IVector<N>& index, const TVector<std::vector<T>, N>& vec) {
  T prod = vec[0][index[0]];
  for(int i = 1; i < N; ++i) prod = prod * vec[i][index[i]];
  return prod;
}

//! Abstract TVector<T, N> from TVector<std::vector<T>, N> by index
/*! e.g. index = { i, j, k, l }
 *  return { vec[0][i], vec[1][j], vec[2][k], vec[3][l] } */
template<typename T, size_t N>
inline TVector<T, N> operator& (const TVector<std::vector<T>, N>& vec, const IVector<N>& index) {
  TVector<T, N> tshape;
  for(int i = 0; i < N; ++i) tshape[i] = vec[i][index[i]];
  return std::move(tshape);
}

//! Abstract TVector<T, N> from TVector<std::vector<T>, N> by index
/*! e.g. index = { i, j, k, l }
 *  return { vec[0][i], vec[1][j], vec[2][k], vec[3][l] } */
template<typename T, size_t N>
inline TVector<T, N> operator& (const IVector<N>& index, const TVector<std::vector<T>, N>& vec) {
  TVector<T, N> tshape;
  for(int i = 0; i < N; ++i) tshape[i] = vec[i][index[i]];
  return std::move(tshape);
}

//####################################################################################################
// Dot product: dot
//####################################################################################################

//! Dot product of two integer vectors
template<size_t N>
inline int dot(const IVector<N>& v1, const IVector<N>& v2) {
  int idot = 0;
  for(int i = 0; i < N; ++i) idot += v1[i]*v2[i];
  return idot;
}

//! Dot product with small overhead, specialized for N = 1
template<>
inline int dot<1>(const IVector<1>& v1, const IVector<1>& v2) {
  return v1[0]*v2[0];
}

//! Dot product with small overhead, specialized for N = 2
template<>
inline int dot<2>(const IVector<2>& v1, const IVector<2>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1];
}

//! Dot product with small overhead, specialized for N = 3
template<>
inline int dot<3>(const IVector<3>& v1, const IVector<3>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

//! Dot product with small overhead, specialized for N = 4
template<>
inline int dot<4>(const IVector<4>& v1, const IVector<4>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
}

//! Dot product with small overhead, specialized for N = 5
template<>
inline int dot<5>(const IVector<5>& v1, const IVector<5>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4];
}

//! Dot product with small overhead, specialized for N = 6
template<>
inline int dot<6>(const IVector<6>& v1, const IVector<6>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5];
}

//! Dot product with small overhead, specialized for N = 7
template<>
inline int dot<7>(const IVector<7>& v1, const IVector<7>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5]+v1[6]*v2[6];
}

//! Dot product with small overhead, specialized for N = 8
template<>
inline int dot<8>(const IVector<8>& v1, const IVector<8>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5]+v1[6]*v2[6]+v1[7]*v2[7];
}

//! Dot product with small overhead, specialized for N = 9
template<>
inline int dot<9>(const IVector<9>& v1, const IVector<9>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5]+v1[6]*v2[6]+v1[7]*v2[7]+v1[8]*v2[8];
}

//! Dot product with small overhead, specialized for N = 10
template<>
inline int dot<10>(const IVector<10>& v1, const IVector<10>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5]+v1[6]*v2[6]+v1[7]*v2[7]+v1[8]*v2[8]+v1[9]*v2[9];
}

//! Dot product with small overhead, specialized for N = 11
template<>
inline int dot<11>(const IVector<11>& v1, const IVector<11>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5]+v1[6]*v2[6]+v1[7]*v2[7]+v1[8]*v2[8]+v1[9]*v2[9]+v1[10]*v2[10];
}

//! Dot product with small overhead, specialized for N = 12
template<>
inline int dot<12>(const IVector<12>& v1, const IVector<12>& v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5]+v1[6]*v2[6]+v1[7]*v2[7]+v1[8]*v2[8]+v1[9]*v2[9]+v1[10]*v2[10]+v1[11]*v2[11];
}

//####################################################################################################
// IVector constructor: shape
//####################################################################################################

//! Convenient IVector constructor for N = 1
inline IVector< 1> shape(int n01) {
  IVector< 1> _shape = { n01 };
  return _shape;
}

//! Convenient IVector constructor for N = 2
inline IVector< 2> shape(int n01, int n02) {
  IVector< 2> _shape = { n01, n02 };
  return _shape;
}

//! Convenient IVector constructor for N = 3
inline IVector< 3> shape(int n01, int n02, int n03) {
  IVector< 3> _shape = { n01, n02, n03 };
  return _shape;
}

//! Convenient IVector constructor for N = 4
inline IVector< 4> shape(int n01, int n02, int n03, int n04) {
  IVector< 4> _shape = { n01, n02, n03, n04 };
  return _shape;
}

//! Convenient IVector constructor for N = 5
inline IVector< 5> shape(int n01, int n02, int n03, int n04, int n05) {
  IVector< 5> _shape = { n01, n02, n03, n04, n05 };
  return _shape;
}

//! Convenient IVector constructor for N = 6
inline IVector< 6> shape(int n01, int n02, int n03, int n04, int n05, int n06) {
  IVector< 6> _shape = { n01, n02, n03, n04, n05, n06 };
  return _shape;
}

//! Convenient IVector constructor for N = 7
inline IVector< 7> shape(int n01, int n02, int n03, int n04, int n05, int n06, int n07) {
  IVector< 7> _shape = { n01, n02, n03, n04, n05, n06, n07 };
  return _shape;
}

//! Convenient IVector constructor for N = 8
inline IVector< 8> shape(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08) {
  IVector< 8> _shape = { n01, n02, n03, n04, n05, n06, n07, n08 };
  return _shape;
}

//! Convenient IVector constructor for N = 9
inline IVector< 9> shape(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09) {
  IVector< 9> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09 };
  return _shape;
}

//! Convenient IVector constructor for N = 10
inline IVector<10> shape(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10) {
  IVector<10> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10 };
  return _shape;
}

//! Convenient IVector constructor for N = 11
inline IVector<11> shape(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11) {
  IVector<11> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11 };
  return _shape;
}

//! Convenient IVector constructor for N = 12
inline IVector<12> shape(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11, int n12) {
  IVector<12> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11, n12 };
  return _shape;
}

//! Convenient IVector constructor to return sequence
/*! e.g. sequence<5>(0, 1) returns { 0, 1, 2, 3, 4 }
 *  \param first starting value of sequence
 *  \param incl  increments of sequence */
template<size_t N>
IVector<N> sequence(int first = 0, int incl = 1) {
  IVector<N> seq;
  // Explicit fixed-size loop might be faster than using std::iota?
  for(int i = 0; i < N; ++i) {
    seq[i] = first;
    first += incl;
  }
  return std::move(seq);
}

//####################################################################################################
// Transpose and Permute vector elements
//####################################################################################################

//! Return transposed vector.
/*! e.g. N = 6, K = 4
 *  [i,j,k,l,m,n] -> [m,n,i,j,k,l]
 *  |1 2 3 4|5 6|    |5 6|1 2 3 4| */
template<typename T, size_t N>
TVector<T, N> transpose(const TVector<T, N>& vec, int K) {
  assert(K >= 0 && K <= N);
  TVector<T, N> vtr;
  for(int i = 0; i < N-K; ++i) vtr[i] = vec[i+K];
  for(int i = N-K; i < N; ++i) vtr[i] = vec[i+K-N];
  return std::move(vtr);
}

//! Return permuted vector
template<typename T, size_t N>
TVector<T, N> permute(const TVector<T, N>& vec, const IVector<N>& pindex) {
  TVector<T, N> vpm;
  for(int i = 0; i < N; ++i) vpm[i] = vec[pindex[i]];
  return std::move(vpm);
}

//####################################################################################################
// Direct product of Dshapes ( aka std::vector<int> ) as operator*
//####################################################################################################

//! Direct product of Dshapes
inline Dshapes operator* (const Dshapes& ds1, const Dshapes& ds2) {
  Dshapes dpr;
  dpr.reserve(ds1.size()*ds2.size());
  for(const int& di : ds1)
    for(const int& dj : ds2) dpr.push_back(di*dj);
  return std::move(dpr);
}

//####################################################################################################
// Fast copying and adding function: fast_copy, fast_add
//####################################################################################################

//! Fast copy for general type
template<typename T>
inline void _fast_copy(size_t n, const T* x, T* y) { std::copy(x, x+n, y); }

//! Fast scale function for general type (this must be specialized).
template<typename T>
inline void _fast_scal(size_t n, const T& alpha, T* x) { }

//! Fast addition for general type (this must be specialized).
template<typename T>
inline void _fast_add (size_t n, const T* x, T* y) { }

//! Fast copy function for std::vector<T>
template<typename T>
void fast_copy(const std::vector<T>& v1, std::vector<T>& v2) {
  size_t n = v1.size(); v2.resize(n);
  // specialized by T
  _fast_copy(n, v1.data(), v2.data());
}

//! Fast scale function for std::vector<T>
template<typename T>
void fast_scal(const T& alpha, std::vector<T>& v) {
  size_t n = v.size();
  // specialized by T
  _fast_scal(n, alpha, v.data());
}

//! Fast copy function for std::vector<T>
template<typename T>
void fast_add(const std::vector<T>& v1, std::vector<T>& v2) {
  size_t n1 = v1.size();
  size_t n2 = v2.size();
  assert(n1 == n2);
  // specialized by T
  _fast_add(n1, v1.data(), v2.data());
}

}; // namespace btas

//####################################################################################################
// Stream printing operator for std::vector<T>
//####################################################################################################

#include <iostream>

//! Printing elements in array as " [ v[0], v[1], v[2], ... ] "
template<typename T, size_t N>
std::ostream& operator<< (std::ostream& ost, const std::array<T, N>& vec) {
  ost << "[ ";
  for(int i = 0; i < N-1; ++i) ost << vec[i] << ", ";
  ost << vec[N-1] << " ]";
  return ost;
}

//! Printing elements in vector as " [ v[0], v[1], v[2], ... ] "
template<typename T>
std::ostream& operator<< (std::ostream& ost, const std::vector<T>& vec) {
  int n = vec.size();
  if(n == 0) {
    ost << "[ ]";
  }
  else {
    ost << "[ ";
    for(int i = 0; i < n-1; ++i) ost << vec[i] << ", ";
    ost << vec[n-1] << " ]";
  }
  return ost;
}

#endif // _BTAS_CXX11_TVECTOR_H
