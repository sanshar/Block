
//####################################################################################################
/*! \file  QSDblas.h 
 *  \brief BLAS wrappers for quantum number-based sparse array contractions */
//####################################################################################################

#ifndef _BTAS_CXX11_QSDBLAS_H
#define _BTAS_CXX11_QSDBLAS_H 1

#include <btas/btas.h>
#include <btas/btas_contract_shape.h>

#include <btas/SPARSE/SDblas.h>

#include <btas/QSPARSE/btas_contract_qshape.h>
#include <btas/QSPARSE/QSDArray.h>

namespace btas {

//====================================================================================================
// BLAS LEVEL 1
//====================================================================================================

//! DCOPY
template<size_t N, class Q = Quantum>
void QSDcopy
(const QSDArray<N, Q>& x, QSDArray<N, Q>& y)
{
  y.resize(x.q(), x.qshape(), x.dshape(), false);
#ifdef _SERIAL
  serial_SDcopy(x, y, false);
#else
  thread_SDcopy(x, y, false);
#endif
}

//! DSCAL
template<size_t N, class Q = Quantum>
void QSDscal
(const double& alpha, QSDArray<N, Q>& x)
{
  SDscal(alpha, x);
}

//! DDOT(U): X^T * Y
template<size_t N, class Q = Quantum>
double QSDdotu
(const QSDArray<N, Q>& x, const QSDArray<N, Q>& y)
{
  if(x.q() != -y.q())
    BTAS_THROW(false, "btas::QSDdotu: total quantum number mismatched");
  if(x.qshape() != -y.qshape())
    BTAS_THROW(false, "btas::QSDdotu: quantum number indices mismatched");
  return serial_SDdot(x, y);
}

//! DDOT(C): X^H * Y
//! Taking hermitian conjugate of x upon contraction (checking quantum number process only differs from QSDdotu)
template<size_t N, class Q = Quantum>
double QSDdotc
(const QSDArray<N, Q>& x, const QSDArray<N, Q>& y)
{
  if(x.q() != y.q())
    BTAS_THROW(false, "btas::QSDdotc: total quantum number mismatched");
  if(x.qshape() != y.qshape())
    BTAS_THROW(false, "btas::QSDdotc: quantum number indices mismatched");
  return serial_SDdot(x, y);
}

//! DAXPY
template<size_t N, class Q = Quantum>
void QSDaxpy
(const double& alpha, const QSDArray<N, Q>& x, QSDArray<N, Q>& y)
{
  if(y.size() > 0) {
    if(x.q() != y.q())
      BTAS_THROW(false, "btas::QSDaxpy: total quantum number mismatched");
    if(x.qshape() != y.qshape())
      BTAS_THROW(false, "btas::QSDaxpy: quantum number indices mismatched");
    if(x.dshape() != y.dshape())
      BTAS_THROW(false, "btas::SDaxpy: shape of y mismatched");
  }
  else {
    y.resize(x.q(), x.qshape(), x.dshape(), false);
  }
#ifdef _SERIAL
  serial_SDaxpy(alpha, x, y);
#else
  thread_SDaxpy(alpha, x, y);
#endif
}

//====================================================================================================
// BLAS LEVEL 2
//====================================================================================================

//! DGEMV
template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDgemv
(const BTAS_TRANSPOSE& TransA,
 const double& alpha, const QSDArray<NA, Q>& a, const QSDArray<NB, Q>& b,
 const double& beta,        QSDArray<NC, Q>& c)
{
  // Checking contraction quantum numbers
  Q q_total;
  TVector<Qshapes<Q>, NC> q_shape;
  gemv_contract_qshape(TransA, a.q(), a.qshape(), b.q(), b.qshape(), q_total, q_shape);
  TVector<Dshapes,    NC> d_shape;
  gemv_contract_dshape(TransA,        a.dshape(),        b.dshape(),          d_shape);
  // Checking shapes of c
  if(c.size() > 0) {
    if(q_total != c.q())
      BTAS_THROW(false, "btas::QSDgemv: total quantum number of c mismatched");
    if(q_shape != c.qshape())
      BTAS_THROW(false, "btas::QSDgemv: quantum number indices of c mismatched");
    if(d_shape != c.dshape())
      BTAS_THROW(false, "btas::QSDgemv: block shape of c mismatched");
    SDscal(beta, c);
  }
  else {
    c.resize(q_total, q_shape, d_shape, false);
  }
  // Calling SDgemv
  if(TransA == NoTrans)
    thread_SDgemv(TransA, alpha, a, b, c);
  else
    thread_SDgemv(TransA, alpha, a.transposed_view(NB), b, c);
}

//! DGER
template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDger
(const double& alpha, const QSDArray<NA, Q>& a, const QSDArray<NB, Q>& b, QSDArray<NC, Q>& c)
{
  // Checking contraction quantum numbers
  Q q_total;
  TVector<Qshapes<Q>, NC> q_shape;
  ger_contract_qshape(a.q(), a.qshape(), b.q(), b.qshape(), q_total, q_shape);
  TVector<Dshapes,    NC> d_shape;
  ger_contract_dshape(       a.dshape(),        b.dshape(),          d_shape);
  // Checking shapes of c
  if(c.size() > 0) {
    if(q_total != c.q())
      BTAS_THROW(false, "btas::QSDger: total quantum number of c mismatched");
    if(q_shape != c.qshape())
      BTAS_THROW(false, "btas::QSDger: quantum number indices of c mismatched");
    if(d_shape != c.dshape())
      BTAS_THROW(false, "btas::QSDger: block shape of c mismatched");
  }
  else {
    c.resize(q_total, q_shape, d_shape, false);
  }
  // Calling SDger
  thread_SDger(alpha, a, b, c);
}

//====================================================================================================
// BLAS LEVEL 3
//====================================================================================================

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDgemm
(const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
 const double& alpha, const QSDArray<NA, Q>& a, const QSDArray<NB, Q>& b,
 const double& beta,        QSDArray<NC, Q>& c)
{
  // Checking contraction quantum numbers
  const size_t K = (NA + NB - NC)/2;
  Q q_total;
  TVector<Qshapes<Q>, NC> q_shape;
  gemm_contract_qshape(TransA, TransB, a.q(), a.qshape(), b.q(), b.qshape(), q_total, q_shape);
  TVector<Dshapes,    NC> d_shape;
  gemm_contract_dshape(TransA, TransB,        a.dshape(),        b.dshape(),          d_shape);
  // Checking shapes of c
  if(c.size() > 0) {
    if(q_total != c.q())
      BTAS_THROW(false, "btas::QSDgemm: total quantum number of c mismatched");
    if(q_shape != c.qshape())
      BTAS_THROW(false, "btas::QSDgemm: quantum number indices of c mismatched");
    if(d_shape != c.dshape())
      BTAS_THROW(false, "btas::QSDgemm: block shape of c mismatched");
    SDscal(beta, c);
  }
  else {
    c.resize(q_total, q_shape, d_shape, false);
  }
  // Calling SDgemm
  if     (TransA == NoTrans && TransB == NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a, b.transposed_view(K), c);

  else if(TransA == NoTrans && TransB != NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a, b, c);

  else if(TransA != NoTrans && TransB == NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a.transposed_view(K), b.transposed_view(K), c);

  else if(TransA != NoTrans && TransB != NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a.transposed_view(K), b, c);
}

//####################################################################################################
// QSD-BLAS with index-based contraction scaling s.t. Clebsch-Gordan coefficient
//####################################################################################################

//====================================================================================================
// BLAS LEVEL 2
//====================================================================================================

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDgemv
(const function<double(const TVector<Q, NA>&,
                       const TVector<Q, NB>&,
                       const TVector<Q, NC>&)>& f_scale,
 const BTAS_TRANSPOSE& TransA,
 const double& alpha, const QSDArray<NA>& a, const QSDArray<NB>& b,
 const double& beta,        QSDArray<NC>& c)
{
  // Checking contraction quantum numbers
  Q q_total;
  TVector<Qshapes<Q>, NC> q_shape;
  gemv_contract_qshape(TransA, a.q(), a.qshape(), b.q(), b.qshape(), q_total, q_shape);
  TVector<Dshapes,    NC> d_shape;
  gemv_contract_dshape(TransA,        a.dshape(),        b.dshape(),          d_shape);
  // Checking shapes of c
  if(c.size() > 0) {
    if(q_total != c.q())
      BTAS_THROW(false, "btas::QSDgemv: total quantum number of c mismatched");
    if(q_shape != c.qshape())
      BTAS_THROW(false, "btas::QSDgemv: quantum number indices of c mismatched");
    if(d_shape != c.dshape())
      BTAS_THROW(false, "btas::QSDgemv: block shape of c mismatched");
    SDscal(beta, c);
  }
  else {
    c.resize(q_total, q_shape, d_shape, false);
  }
  // Calling SDgemv with index-based scaling
  if(TransA == NoTrans) {
    function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)>
    q_scale = bind(f_indexbase_scale<double, NA, NB, NC, Q>, f_scale, a.qshape(), _1, b.qshape(), _2, c.qshape(), _3);
    thread_SDgemv(q_scale, TransA, alpha, a, b, c);
  }
  else {
    function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)>
    // FIXME: Is this correct? Transposition may introduce change something more?
    q_scale = bind(f_indexbase_scale<double, NA, NB, NC, Q>, f_scale, transpose(a.qshape(), NB), _1, b.qshape(), _2, c.qshape(), _3);
    thread_SDgemv(q_scale, TransA, alpha, a.transposed_view(NB), b, c);
  }
}

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDger
(const function<double(const TVector<Q, NA>&,
                       const TVector<Q, NB>&,
                       const TVector<Q, NC>&)>& f_scale,
 const double& alpha, const QSDArray<NA>& a, const QSDArray<NB>& b, QSDArray<NC>& c)
{
  // Checking contraction quantum numbers
  Q q_total;
  TVector<Qshapes<Q>, NC> q_shape;
  ger_contract_qshape(a.q(), a.qshape(), b.q(), b.qshape(), q_total, q_shape);
  TVector<Dshapes,    NC> d_shape;
  ger_contract_dshape(       a.dshape(),        b.dshape(),          d_shape);
  // Checking shapes of c
  if(c.size() > 0) {
    if(q_total != c.q())
      BTAS_THROW(false, "btas::QSDger: total quantum number of c mismatched");
    if(q_shape != c.qshape())
      BTAS_THROW(false, "btas::QSDger: quantum number indices of c mismatched");
    if(d_shape != c.dshape())
      BTAS_THROW(false, "btas::QSDger: block shape of c mismatched");
  }
  else {
    c.resize(q_total, q_shape, d_shape, false);
  }
  // Calling SDger with index-based scaling
  function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)>
  q_scale = bind(f_indexbase_scale<double, NA, NB, NC, Q>, f_scale, a.qshape(), _1, b.qshape(), _2, c.qshape(), _3);
  thread_SDger(q_scale, alpha, a, b, c);
}

//====================================================================================================
// BLAS LEVEL 3
//====================================================================================================

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDgemm
(const function<double(const TVector<Q, NA>&,
                       const TVector<Q, NB>&,
                       const TVector<Q, NC>&)>& f_scale,
 const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
 const double& alpha, const QSDArray<NA>& a, const QSDArray<NB>& b,
 const double& beta,        QSDArray<NC>& c)
{
  // Checking contraction quantum numbers
  const size_t K = (NA + NB - NC)/2;
  Q q_total;
  TVector<Qshapes<Q>, NC> q_shape;
  gemm_contract_qshape(TransA, TransB, a.q(), a.qshape(), b.q(), b.qshape(), q_total, q_shape);
  TVector<Dshapes,    NC> d_shape;
  gemm_contract_dshape(TransA, TransB,        a.dshape(),        b.dshape(),          d_shape);
  // Checking shapes of c
  if(c.size() > 0) {
    if(q_total != c.q())
      BTAS_THROW(false, "btas::QSDgemm: total quantum number of c mismatched");
    if(q_shape != c.qshape())
      BTAS_THROW(false, "btas::QSDgemm: quantum number indices of c mismatched");
    if(d_shape != c.dshape())
      BTAS_THROW(false, "btas::QSDgemm: block shape of c mismatched");
    SDscal(beta, c);
  }
  else {
    c.resize(q_total, q_shape, d_shape, false);
  }
  // Calling SDgemm with index-based scaling
    // FIXME: Are those correct? Transposition may introduce change something more?
  if     (TransA == NoTrans && TransB == NoTrans) {
    function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)>
    q_scale = bind(f_indexbase_scale<double, NA, NB, NC, Q>,
                   f_scale, a.qshape(), _1, transpose(b.qshape(), K), _2, c.qshape(), _3);
    thread_SDgemm(q_scale, TransA, TransB, alpha, a, b.transposed_view(K), c);
  }
  else if(TransA == NoTrans && TransB != NoTrans) {
    function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)>
    q_scale = bind(f_indexbase_scale<double, NA, NB, NC, Q>,
                   f_scale, a.qshape(), _1, b.qshape(), _2, c.qshape(), _3);
    thread_SDgemm(q_scale, TransA, TransB, alpha, a, b, c);
  }
  else if(TransA != NoTrans && TransB == NoTrans) {
    function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)>
    q_scale = bind(f_indexbase_scale<double, NA, NB, NC, Q>,
                   f_scale, transpose(a.qshape(), K), _1, transpose(b.qshape(), K), _2, c.qshape(), _3);
    thread_SDgemm(q_scale, TransA, TransB, alpha, a.transposed_view(K), b.transposed_view(K), c);
  }
  else if(TransA != NoTrans && TransB != NoTrans) {
    function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)>
    q_scale = bind(f_indexbase_scale<double, NA, NB, NC, Q>,
                   f_scale, transpose(a.qshape(), K), _1, b.qshape(), _2, c.qshape(), _3);
    thread_SDgemm(q_scale, TransA, TransB, alpha, a.transposed_view(K), b, c);
  }
}

//====================================================================================================
// Non-BLAS utility functions
//====================================================================================================

//! Normalization
template<size_t N, class Q = Quantum>
void QSDnormalize(QSDArray<N, Q>& x) {
  double norm = QSDdotc(x, x);
  QSDscal(1.0/sqrt(norm), x);
}

//! Orthogonalization
template<size_t N, class Q = Quantum>
void QSDorthogonalize(const QSDArray<N, Q>& x, QSDArray<N, Q>& y) {
  double ovlp = QSDdotc(x, y);
  QSDaxpy(-ovlp, x, y);
}

}; // namespace btas

#endif // _BTAS_CXX11_QSDBLAS_H
