//
/*! \file Dlapack.h
 *  \brief LAPACK wrappers for double precision real array
 * 
 *  Implemented in terms of Row-Major array as C/C++ style
 *  Since they contain transposition of array before and after calling LAPACK subroutines,
 *  performance is not good as FORTRAN
 *
 *  Or, if you work with Col-Major array, use clapack_xxx subroutines directly (see lapack/clapack.h)
 */

#ifndef _BTAS_CXX11_DLAPACK_H
#define _BTAS_CXX11_DLAPACK_H 1

#include <algorithm>
#include <numeric>

#include <lapack/clapack.h>

#include <btas/btas.h>
#include <btas/TVector.h>

#include <btas/DENSE/DArray.h>

namespace btas {

//! LU decomposition for general array
/*! i.e. compute fragmentation s.t. A = P * L * U
 *  \param a is an array viewed as a general rectangular matrix A
 *  \param l is an array viewed as a lower triangular matrix L
 *  \param u is an array viewed as a upper triangular matrix U
 *  \param ipiv contains pivot info. corresponding permutation matrix P */
template<size_t M, size_t N>
void Dgetrf(const DArray<M+N-2>& a, DArray<M>& l, DArray<N>& u, std::vector<int>& ipiv) {
  // calc. shapes
  const IVector<M+N-2>& a_shape = a.shape();
  int nrows = std::accumulate(a_shape.begin(),     a_shape.begin()+M-1, 1, std::multiplies<int>());
  int ncols = std::accumulate(a_shape.begin()+M-1, a_shape.end(),       1, std::multiplies<int>());
  int nsize = std::min(nrows, ncols);
  ipiv.resize(nsize);
  // call dgetrf
  DArray<M+N-2> a_tr;
  Dpermute(a, transpose(sequence<M+N-2>(), M-1), a_tr);
  if(clapack_dgetrf(ncols, nrows, a_tr.data(), nrows, ipiv.data()) < 0)
    BTAS_THROW(false, "btas::Dgetrf terminated abnormally");
  // abstract matrix elements
  DArray<2> l_matrix(nrows, nsize); l_matrix = 0.0;
  DArray<2> u_matrix(nsize, ncols); u_matrix = 0.0;
  double* p = a_tr.data();
  for(int j = 0; j < ncols; ++j) {
    for(int i = 0; i < nrows; ++i, ++p) {
      if(i == j)
        l_matrix(i, j) = 1.0;
      if(i >  j)
        l_matrix(i, j) = *p;
      else
        u_matrix(i, j) = *p;
    }
  }
  // reshape 'l' from matrix to tensor
  IVector<M> l_shape;
  for(int i = 0; i < M-1; ++i) l_shape[i] = a_shape[i];
  l_shape[M-1] = nsize;
  l.resize(l_shape);
  DcopyFlatten(l_matrix, l);
  // reshape 'u' from matrix to tensor
  IVector<N> u_shape;
  u_shape[0]   = nsize;
  for(int i = 1; i < N;   ++i) u_shape[i] = a_shape[i+M-2];
  u.resize(u_shape);
  DcopyFlatten(u_matrix, u);
}

//! Solve linear equation A*x = b
/*! \param a is [(M-N+1) x (N-1)] array viewed as a square matrix A
 *  \param b is [(N-1) x (1)] array viewed as a set of input vectors B
 *  \param x is a set of solutions */
template<size_t M, size_t N>
void Dgesv(const DArray<M>& a, const DArray<N>& b, DArray<N>& x) {
  const size_t C = N-1;
  const size_t R = M-C; assert(R > 0);
  if(a.size() == 0)
    BTAS_THROW(false, "btas::Dgesv: array data of 'a' not found");
  if(b.size() == 0)
    BTAS_THROW(false, "btas::Dgesv: array data of 'b' not found");
  // check array shapes of 'a' and 'b'
  const IVector<M>& a_shape = a.shape();
        IVector<N>  b_shape;
  int nrows, ncols, nsols;
  if(C > 0) {
    nrows = std::accumulate(a_shape.begin(),     a_shape.begin()+R,   1, std::multiplies<int>());
    ncols = std::accumulate(a_shape.begin()+R,   a_shape.end(),       1, std::multiplies<int>());
    for(int i = 0; i < C; ++i) b_shape[i] = a_shape[R+i];
    if(!std::equal(b_shape.begin(), b_shape.begin()+C, b.shape().begin()))
      BTAS_THROW(false, "btas::Dgesv: array shape of 'b' mismatched");
    nsols = b.shape(C);
  }
  else {
    nrows = std::accumulate(a_shape.begin(),     a_shape.begin()+R-1, 1, std::multiplies<int>());
    ncols = std::accumulate(a_shape.begin()+R-1, a_shape.end(),       1, std::multiplies<int>());
    if(a.shape(M-1) != b.shape(0))
      BTAS_THROW(false, "btas::Dgesv: array shape of 'b' mismatched");
    nsols = 1;
  }
  // nrows: number of equations
  // ncols: number of variables
  assert(nrows == ncols);
  // transpose 'a' and 'b' (row-major (C/C++) -> col-major (FORTRAN))
  DArray<M> a_tr;
  DArray<N> x_tr;
  if(C > 0) {
    Dpermute(a, transpose(sequence<M>(), R),   a_tr);
    Dpermute(b, transpose(sequence<N>(), C),   x_tr);
  }
  else {
    Dpermute(a, transpose(sequence<M>(), R-1), a_tr);
    Dcopy(b, x_tr);
  }
  std::vector<int> ipiv(nrows);
  // call dgesv
  if(clapack_dgesv(nrows, nsols, a_tr.data(), nrows, ipiv.data(), x_tr.data(), ncols) < 0)
    BTAS_THROW(false, "btas::Dgesv terminated abnormally");
  // transpose back for solution array 'x'
  if(C > 0)
    Dpermute(x_tr, transpose(sequence<N>(), 1), x);
  else
    Dcopy(x_tr, x);
}

//! Solve eigenvalue problem for real symmetric matrix (SEP)
/*! i.e. compute A = Z * D * Z^T
 *  \param a is [(N-1) x (N-1)] symmetric array A
 *  \param d is diagonal matrix D containing eigenvalues
 *  \param z is [(N-1) x (1)] unitary array Z containing eigenvectors */
template<size_t N>
void Dsyev(const DArray<2*N-2>& a, DArray<1>& d, DArray<N>& z) {
  const size_t M = 2*N-2;
  const size_t C = N-1;
  const IVector<M>& a_shape = a.shape();
  // check array 'a'
  if(a.size() == 0)
    BTAS_THROW(false, "btas::Dsyev: array data not found");
  if(!std::equal(a_shape.begin(), a_shape.begin()+C, a_shape.begin()+C))
    BTAS_THROW(false, "btas::Dsyev: array shape of 'a' must be symmetric");
  // calc. column shape and its total size
  IVector<C>  c_shape;
  for(int i = 0; i < C; ++i) c_shape[i] = a_shape[i];
  int ncols = std::accumulate(c_shape.begin(), c_shape.end(), 1, std::multiplies<int>());
  // z_shape is transposed shape of eigenvectors
  IVector<N> z_shape;
  z_shape[0] = ncols;
  for(int i = 1; i < N; ++i) z_shape[i] = c_shape[i-1];
  // assumed 'a' is symmetric
  DArray<N> z_tr(z_shape);
  DcopyFlatten(a, z_tr);
  // resize eigenvalue vector
  d.resize(ncols);
  // call dsyev
  if(clapack_dsyev(ClapackCalcVector, ClapackUpper, ncols, z_tr.data(), ncols, d.data()) < 0)
    BTAS_THROW(false, "btas::Dsyev terminated abnormally");
  // transpose eigenvectors to [(N-1) x (1)] array
  Dpermute(z_tr, transpose(sequence<N>(), 1), z);
}

//! Solve eigenvalue problem for real non-symmetric square matrix (NSGEP)
/*! [(N-1) x (N-1)] square array 'a' is decomposed into [(N-1) x (1)] left/right eigenvectors
 *  and corresponding complex diagonal matrix
 *
 *  i.e. [ A * vr(j) = w(j) * vr(j) ] and [ vl(j)^T * A = w(j) * vl(j)^T ], where [ 'w' = 'wr' + i * 'wi' ]
 *
 *  \param a  array to be decomposed
 *  \param wr real part of eigenvalues
 *  \param wi imaginary part of eigenvalues
 *  \param vl left eigenvectors
 *  \param vr right eigenvectors
 */
template<size_t N>
void Dgeev(const DArray<2*N-2>& a, DArray<1>& wr, DArray<1>& wi, DArray<N>& vl, DArray<N>& vr) {
  const size_t M = 2*N-2;
  const size_t C = N-1;
  const IVector<M>& a_shape = a.shape();
  // check array 'a'
  if(a.size() == 0)
    BTAS_THROW(false, "btas::Dgeev: array data not found");
  if(!std::equal(a_shape.begin(), a_shape.begin()+C, a_shape.begin()+C))
    BTAS_THROW(false, "btas::Dgeev: array shape of 'a' must be symmetric");
  DArray<M> a_tr;
  Dpermute(a, transpose(sequence<M>(), C), a_tr);
  // calc. column shape and its total size
  IVector<C> c_shape;
  for(int i = 0; i < C; ++i) c_shape[i] = a_shape[i];
  int ncols = std::accumulate(c_shape.begin(), c_shape.end(), 1, std::multiplies<int>());
  // v_shape is transposed shape of eigenvectors
  IVector<N> v_shape;
  v_shape[0] = ncols;
  for(int i = 1; i < N; ++i) v_shape[i] = c_shape[i-1];
  DArray<N> vl_tr(v_shape);
  DArray<N> vr_tr(v_shape);
  // resize eigenvalue vector
  wr.resize(ncols);
  wi.resize(ncols);
  // call dgeev
  if(clapack_dgeev(ClapackCalcVector, ClapackCalcVector,
                   ncols, a_tr.data(), ncols, wr.data(), wi.data(), vl_tr.data(), ncols, vr_tr.data(), ncols) < 0)
    BTAS_THROW(false, "btas::Dgeev terminated abnormally");
  // transpose eigenvectors to [(N-1) x (1)] array
  Dpermute(vl_tr, transpose(sequence<N>(), 1), vl);
  Dpermute(vr_tr, transpose(sequence<N>(), 1), vr);
}

//! Solve generalized eigenvalue problem for real symmetric matrix
/*! i.e. [ A * z(j) = d(j) * B * z(j) ] */
template<size_t N>
void Dsygv(const DArray<2*N-2>& a, const DArray<2*N-2>& b, DArray<1>& d, DArray<N>& z) {
  const size_t M = 2*N-2;
  const size_t C = N-1;
  const IVector<M>& a_shape = a.shape();
  const IVector<M>& b_shape = b.shape();
  if(a_shape != b_shape)
    BTAS_THROW(false, "btas::Dsygv: shapes of 'a' and 'b' mismatched");
  // check array 'a'
  if(a.size() == 0)
    BTAS_THROW(false, "btas::Dsygv: array data of 'a' not found");
  if(!std::equal(a_shape.begin(), a_shape.begin()+C, a_shape.begin()+C))
    BTAS_THROW(false, "btas::Dsygv: array shape of 'a' must be symmetric");
  // check array 'b'
  if(b.size() == 0)
    BTAS_THROW(false, "btas::Dsygv: array data of 'b' not found");
  if(!std::equal(b_shape.begin(), b_shape.begin()+C, b_shape.begin()+C))
    BTAS_THROW(false, "btas::Dsygv: array shape of 'b' must be symmetric");
  // calc. column shape and its total size
  IVector<C> c_shape;
  for(int i = 0; i < C; ++i) c_shape[i] = a_shape[i];
  int ncols = std::accumulate(c_shape.begin(), c_shape.end(), 1, std::multiplies<int>());
  // z_shape is transposed shape of eigenvectors
  IVector<N> z_shape;
  z_shape[0] = ncols;
  for(int i = 1; i < N; ++i) z_shape[i] = c_shape[i-1];
  // assumed 'a' is symmetric
  DArray<N> z_tr(z_shape);
  DcopyFlatten(a, z_tr);
  // take copy of 'b'
  DArray<M> b_tr(b);
  // resize eigenvalue vector
  d.resize(ncols);
  // call dsygv
  if(clapack_dsygv(1, ClapackCalcVector, ClapackUpper, ncols, z_tr.data(), ncols, b_tr.data(), ncols, d.data()) < 0)
    BTAS_THROW(false, "btas::Dsygv terminated abnormally");
  // transpose eigenvectors to [(N-1) x (1)] array
  Dpermute(z_tr, transpose(sequence<N>(), 1), z);
}

//! Solve generalized eigenvalue problem for non-symmetric real matrix pair
/*! i.e. [ A * vr(j) = d(j) * B * vr(j) ] and [ vl(j)^T * A = d(j) * vl(j)^T * B ]
 *  where [ d(j) = (alphar(j) + i * alphai(j)) / beta(j) ],
 *  and [ d(j+1) = (alphar(j) - i * alphai(j)) / beta(j) ] if alphai(j) is not zero
 */
template<size_t N>
void Dggev(const DArray<2*N-2>& a, const DArray<2*N-2>& b, DArray<1>& alphar, DArray<1>& alphai, DArray<1>& beta, DArray<N>& vl, DArray<N>& vr) {
  const size_t M = 2*N-2;
  const size_t C = N-1;
  const IVector<M>& a_shape = a.shape();
  const IVector<M>& b_shape = b.shape();
  if(a_shape != b_shape)
    BTAS_THROW(false, "btas::Dggev: shapes of 'a' and 'b' mismatched");
  // check array 'a'
  if(a.size() == 0)
    BTAS_THROW(false, "btas::Dggev: array data of 'a' not found");
  if(!std::equal(a_shape.begin(), a_shape.begin()+C, a_shape.begin()+C))
    BTAS_THROW(false, "btas::Dggev: array shape of 'a' must be symmetric");
  // check array 'b'
  if(b.size() == 0)
    BTAS_THROW(false, "btas::Dggev: array data of 'b' not found");
  if(!std::equal(b_shape.begin(), b_shape.begin()+C, b_shape.begin()+C))
    BTAS_THROW(false, "btas::Dggev: array shape of 'b' must be symmetric");
  // transpose 'a' and 'b'
  DArray<M> a_tr;
  Dpermute(a, transpose(sequence<M>(), C), a_tr);
  DArray<M> b_tr;
  Dpermute(b, transpose(sequence<M>(), C), b_tr);
  // calc. column shape and its total size
  IVector<C> c_shape;
  for(int i = 0; i < C; ++i) c_shape[i] = a_shape[i];
  int ncols = std::accumulate(c_shape.begin(), c_shape.end(), 1, std::multiplies<int>());
  // v_shape is transposed shape of eigenvectors
  IVector<N> v_shape;
  v_shape[0] = ncols;
  for(int i = 1; i < N; ++i) v_shape[i] = c_shape[i-1];
  DArray<N> vl_tr(v_shape);
  DArray<N> vr_tr(v_shape);
  // resize eigenvalue vector
  alphar.resize(ncols);
  alphai.resize(ncols);
  beta.resize(ncols);
  // call dggev
  if(clapack_dggev(ClapackCalcVector, ClapackCalcVector, ncols, a_tr.data(), ncols, b_tr.data(), ncols,
                   alphar.data(), alphai.data(), beta.data(), vl_tr.data(), ncols, vr_tr.data(), ncols) < 0)
    BTAS_THROW(false, "btas::Dggev terminated abnormally");
  // transpose eigenvectors to [(N-1) x (1)] array
  Dpermute(vl_tr, transpose(sequence<N>(), 1), vl);
  Dpermute(vr_tr, transpose(sequence<N>(), 1), vr);
}

//! Solve singular value decomposition
template<size_t NA, size_t NU>
void Dgesvd(const DArray<NA>& a, DArray<1>& s, DArray<NU>& u, DArray<NA-NU+2>& vt, bool calc_full_u = false, bool calc_full_vt = false) {
  const size_t K  = NA - NU + 1;
  const size_t NV = K + 1;
  // check array 'a'
  if(a.size() == 0)
    BTAS_THROW(false, "btas::Dgesvd: array data not found");
  // calc. sizes
  const IVector<NA>& a_shape = a.shape();
  int nrows = std::accumulate(a_shape.begin(),     a_shape.begin()+NA-K, 1, std::multiplies<int>());
  int ncols = std::accumulate(a_shape.begin()+NA-K, a_shape.end(),       1, std::multiplies<int>());
  int nsval = std::min(nrows, ncols);
  // calc. shape of 'u'
  int ucols = calc_full_u ? nrows : nsval;
//int ucols = nsval;
  IVector<NU> u_shape;
  for(int i = 0; i < NA-K; ++i) u_shape[i] = a_shape[i];
  u_shape[NA-K] = ucols;
  // calc. shape of 'vt'
  int vrows = calc_full_vt ? ncols : nsval;
//int vrows = nsval;
  IVector<NV> vt_shape;
  vt_shape[0] = vrows;
  for(int i = 0; i < K; ++i) vt_shape[i+1] = a_shape[i+NA-K];
  // resize arrays
  s.resize(nsval);
  u.resize(u_shape);
  vt.resize(vt_shape);
  // take copy of 'a' since it will be destroyed upon calling dgesvd
  DArray<NA> a_cp(a);
  // call dgesvd as [ A^T = V * S * U^T ] instead of [ A = U * S * V^T ]
  if(clapack_dgesvd(ClapackCalcThinVector, ClapackCalcThinVector,
                    ncols, nrows, a_cp.data(), ncols, s.data(), vt.data(), ncols, u.data(), ucols) < 0)
    BTAS_THROW(false, "btas::Dgesvd terminated abnormally");
}

}; // namespace btas

#endif // _BTAS_CXX11_DLAPACK_H
