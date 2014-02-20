//
/*! \file Dblas.h
 *  \brief BLAS wrappers for double precision real array
 */

#ifndef _BTAS_CXX11_DBLAS_H
#define _BTAS_CXX11_DBLAS_H 1

#include <algorithm>
#include <numeric>

#include <btas/btas.h>
#include <btas/btas_contract_shape.h>
#include <btas/blas_cxx_interface.h>

#include <btas/DENSE/DArray.h>

namespace btas {

//####################################################################################################
// BLAS LEVEL 1
//####################################################################################################

//! DCOPY: y := x
template<size_t N>
void Dcopy(const DArray<N>& x, DArray<N>& y) {
  if(x.size() == 0) {
    y.clear();
  }
  else {
    y.resize(x.shape());
    cblas_dcopy(x.size(), x.data(), 1, y.data(), 1);
  }
}

//! DCOPY as flattened vectors: y := x
/*! x and y must have the same size
 *  ranks of x and y can be varied */
template<size_t NX, size_t NY>
void DcopyFlatten(const DArray<NX>& x, DArray<NY>& y) {
  if(x.size() != y.size())
    BTAS_THROW(false, "btas::DcopyFlatten: inconsistent data size");
  cblas_dcopy(x.size(), x.data(), 1, y.data(), 1);
}

//! reshape rank of array
template<size_t NX, size_t NY>
void Dreshape(const DArray<NX>& x, const IVector<NY>& y_shape, DArray<NY>& y) {
  y.resize(y_shape);
  DcopyFlatten(x, y);
}

//! DSCAL: x := alpha * x
template<size_t N>
void Dscal(const double& alpha, DArray<N>& x) {
  if(x.size() == 0) return;
  cblas_dscal(x.size(), alpha, x.data(), 1);
}

//! DAXPY: y := alpha * x + y
template<size_t N>
void Daxpy(const double& alpha, const DArray<N>& x, DArray<N>& y) {
  if(x.size() == 0)
    BTAS_THROW(false, "btas::Daxpy: array data not found");
  if(y.size() > 0) {
    if(x.shape() != y.shape())
      BTAS_THROW(false, "btas::Daxpy: array shape mismatched");
  }
  else {
    y.resize(x.shape());
    y = 0.0;
  }
  cblas_daxpy(x.size(), alpha, x.data(), 1, y.data(), 1);
}

//! DDOT: inner product of x and y
template<size_t N>
double Ddot(const DArray<N>& x, const DArray<N>& y) {
  if(x.shape() != y.shape())
    BTAS_THROW(false, "btas::Ddot: array shape mismatched");
  return cblas_ddot(x.size(), x.data(), 1, y.data(), 1);
}

//! DNRM2: returns euclidean norm of x
template<size_t N>
double Dnrm2(const DArray<N>& x) {
  return cblas_dnrm2(x.size(), x.data(), 1);
}

//####################################################################################################
// BLAS LEVEL 2
//####################################################################################################

//! DGEMV: Matrix-vector multiplication, c := alphe * a * b + beta * c
template<size_t NA, size_t NB, size_t NC>
void Dgemv(const BTAS_TRANSPOSE& TransA,
           const double& alpha, const DArray<NA>& a, const DArray<NB>& b, const double& beta, DArray<NC>& c) {
  if(a.size() == 0 || b.size() == 0)
    BTAS_THROW(false, "btas::Dgemv: array data not found");

  IVector<NC> c_shape;
  gemv_contract_shape(TransA, a.shape(), b.shape(), c_shape);
  // check and resize c
  if(c.size() > 0) {
    if(c_shape != c.shape())
      BTAS_THROW(false, "btas::Dgemv: array shape of c mismatched");
  }
  else {
    c.resize(c_shape);
    c = 0.0;
  }
  // calling cblas_dgemv
  int arows = std::accumulate(c_shape.begin(), c_shape.end(), 1, std::multiplies<int>());
  int acols = a.size() / arows;
  if(TransA != NoTrans) std::swap(arows, acols);
  cblas_dgemv(RowMajor, TransA, arows, acols, alpha, a.data(), acols, b.data(), 1, beta, c.data(), 1);
}

//! DGER: Rank-update / direct product, c := alpha * ( a x b )
template<size_t NA, size_t NB, size_t NC>
void Dger(const double& alpha, const DArray<NA>& a, const DArray<NB>& b, DArray<NC>& c) {
  if(a.size() == 0 || b.size() == 0)
    BTAS_THROW(false, "btas::Dger: array data not found");

  IVector<NC>  c_shape;
  ger_contract_shape(a.shape(), b.shape(), c_shape);
  // check and resize c
  if(c.size() > 0) {
    if(c_shape != c.shape())
      BTAS_THROW(false, "btas::Dger: array shape of c mismatched");
  }
  else {
    c.resize(c_shape);
    c = 0.0;
  }
  // calling cblas_dger
  int crows = a.size();
  int ccols = b.size();
  cblas_dger(RowMajor, crows, ccols, alpha, a.data(), 1, b.data(), 1, c.data(), ccols);
}

//####################################################################################################
// BLAS LEVEL 3
//####################################################################################################

//! DGEMM: Matrix-matrix multiplication, c := alpha * a * b + beta * c
/*! \param alpha scalar number
 *  \param beta  scalar number
 *  \param a     array to be contracted, regarded as matrix
 *  \param b     array to be contracted, regarded as matrix
 *  \param c     array to be returned,   regarded as matrix */
template<size_t NA, size_t NB, size_t NC>
void Dgemm(const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
           const double& alpha, const DArray<NA>& a, const DArray<NB>& b, const double& beta, DArray<NC>& c) {
  const size_t K = (NA + NB - NC)/2; /*! ranks to be contracted */
  if(a.size() == 0 || b.size() == 0)
    BTAS_THROW(false, "btas::Dgemm: array data not found");

  IVector<K> contracts;
  IVector<NC> c_shape;
  gemm_contract_shape(TransA, TransB, a.shape(), b.shape(), contracts, c_shape);
  // check and resize c
  if(c.size() > 0) {
    if(c_shape != c.shape())
      BTAS_THROW(false, "btas::Dgemm: array shape of c mismatched");
  }
  else {
    c.resize(c_shape);
    c = 0.0;
  }
  // calling cblas_dgemm
  int arows = std::accumulate(c_shape.begin(), c_shape.begin()+NA-K, 1, std::multiplies<int>());
  int acols = std::accumulate(contracts.begin(), contracts.end(), 1, std::multiplies<int>());
  int bcols = std::accumulate(c_shape.begin()+NA-K, c_shape.end(), 1, std::multiplies<int>());
  int lda = acols; if(TransA != NoTrans) lda = arows;
  int ldb = bcols; if(TransB != NoTrans) ldb = acols;
  cblas_dgemm(RowMajor, TransA, TransB, arows, bcols, acols, alpha, a.data(), lda, b.data(), ldb, beta, c.data(), bcols);
}

//! Non-BLAS function: a := a(general matrix) * b(diagonal matrix)
/*! NB <= NA */
template<size_t NA, size_t NB>
void Ddimd(DArray<NA>& a, const DArray<NB>& b) {
  const IVector<NA>& a_shape = a.shape();
        IVector<NB>  b_shape;
  for(int i = 0; i < NB; ++i) b_shape[i] = a_shape[i+NA-NB];
  if(!std::equal(b_shape.begin(), b_shape.end(), a_shape.begin()+NA-NB))
    BTAS_THROW(false, "Ddimd: array shape mismatched");

  int nrows = std::accumulate(a_shape.begin(), a_shape.begin()+NA-NB, 1, std::multiplies<int>());
  int ncols = b.size();
  double* pa = a.data();
  for(int i = 0; i < nrows; ++i) {
    const double* pb = b.data();
    for(int j = 0; j < ncols; ++j, ++pa, ++pb)
      (*pa) *= (*pb);
  }
}

//! Non-BLAS function: b := a(diagonal matrix) * b(general matrix)
/*! NA <= NB */
template<size_t NA, size_t NB>
void Ddidm(const DArray<NA>& a, DArray<NB>& b) {
        IVector<NA>  a_shape;
  const IVector<NB>& b_shape = b.shape();
  for(int i = 0; i < NA; ++i) a_shape[i] = b_shape[i];
  if(!std::equal(a_shape.begin(), a_shape.end(), b_shape.begin()))
    BTAS_THROW(false, "Ddidm: array shape mismatched");

  int nrows = a.size();
  int ncols = std::accumulate(b_shape.begin()+NA, b_shape.end(), 1, std::multiplies<int>());
  const double* pa = a.data();
        double* pb = b.data();
  for(int i = 0; i < nrows; ++i, ++pa, pb += ncols)
    cblas_dscal(ncols, *pa, pb, 1);
}

//! Wrapper function for DBLAS
/*! When calling as vector-matrix multiplication (i.e. call DGEMV with Trans),
 *  tranposed array is returned: c = (a * b)^T = b^T * a^T */
template<size_t NA, size_t NB, size_t NC>
void DblasWrapper(const double& alpha, const DArray<NA>& a, const DArray<NB>& b, const double& beta, DArray<NC>& c) {
  const size_t K = (NA + NB - NC)/2;
  if(NA == K) {
    // FIXME: this may be wrong or give undesired result
    //        array c needs to be permuted before and after calling Dgemv
    //        otherwise, user has to take care of it
    Dgemv(Trans,   alpha, b, a, beta, c);
  }
  else if(NB == K) {
    Dgemv(NoTrans, alpha, a, b, beta, c);
  }
  else {
    Dgemm(NoTrans, NoTrans, alpha, a, b, beta, c);
  }
}

//! Normalization
template<size_t N>
void Dnormalize(DArray<N>& x) {
  double nrm2 = Dnrm2(x);
  Dscal(1.0/nrm2, x);
}

//! Orthogonalization
template<size_t N>
void Dorthogonalize(const DArray<N>& x, DArray<N>& y) {
  double ovlp = Ddot(x, y);
  Daxpy(-ovlp, x, y);
}

}; // namespace btas

#endif // _BTAS_CXX11_DBLAS_H
