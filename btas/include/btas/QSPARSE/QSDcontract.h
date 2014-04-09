
//
/*! \file  QSDcontract.h
 *  \brief Convenient contraction functions for QSDArray
 */

#ifndef _BTAS_CXX11_QSDCONTRACT_H
#define _BTAS_CXX11_QSDCONTRACT_H 1

#include <btas/btas_contract_shape.h>

#include <btas/QSPARSE/QSDArray.h>
#include <btas/QSPARSE/QSDblas.h>
#include <btas/QSPARSE/QSDpermute.h>

namespace btas {

template<size_t NA, size_t NB, size_t NC, size_t K, class Q = Quantum>
void QSDcontract
(const double& alpha, const QSDArray<NA, Q>& a, const IVector<K>& a_contract,
                      const QSDArray<NB, Q>& b, const IVector<K>& b_contract,
 const double& beta,        QSDArray<NC, Q>& c)
{
  IVector<NA> a_permute;
  IVector<NB> b_permute;
  unsigned int jobs = get_contract_jobs(a.shape(), a_contract, a_permute,
                                        b.shape(), b_contract, b_permute);

  QSDArray<NA, Q> a_ref;
  if(jobs & JOBMASK_A_PMUTE) QSDpermute(a, a_permute, a_ref);
  else                       a_ref.reference(a);

  BTAS_TRANSPOSE TransA;
  if(jobs & JOBMASK_A_TRANS) TransA = Trans;
  else                       TransA = NoTrans;

  QSDArray<NB, Q> b_ref;
  if(jobs & JOBMASK_B_PMUTE) QSDpermute(b, b_permute, b_ref);
  else                       b_ref.reference(b);

  BTAS_TRANSPOSE TransB;
  if(jobs & JOBMASK_B_TRANS) TransB = Trans;
  else                       TransB = NoTrans;

  switch(jobs & JOBMASK_BLAS_TYPE) {
    case(0):
      QSDgemv(TransA, alpha, a_ref, b_ref, beta, c);
      break;
    case(1):
      QSDgemv(TransB, alpha, b_ref, a_ref, beta, c);
      break;
    case(2):
      QSDgemm(TransA, TransB, alpha, a_ref, b_ref, beta, c);
      break;
    default:
      BTAS_THROW(false, "btas::QSDcontract unknown blas job type returned");
  }
}

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDindexed_contract
(const double& alpha, const QSDArray<NA, Q>& a, const IVector<NA>& a_symbols,
                      const QSDArray<NB, Q>& b, const IVector<NB>& b_symbols,
 const double& beta,        QSDArray<NC, Q>& c, const IVector<NC>& c_symbols)
{
  const size_t K = (NA + NB - NC)/2;
  IVector<K> a_contract;
  IVector<K> b_contract;
  IVector<NC> axb_symbols;
  indexed_contract_shape(a_symbols, a_contract, b_symbols, b_contract, axb_symbols);

  if(c_symbols == axb_symbols) {
    QSDcontract(alpha, a, a_contract, b, b_contract, beta, c);
  }
  else {
    QSDArray<NC, Q> axb;
    if(c.size() > 0)
      QSDindexed_permute(c, c_symbols, axb, axb_symbols);

    QSDcontract(alpha, a, a_contract, b, b_contract, beta, axb);
    QSDindexed_permute(axb, axb_symbols, c, c_symbols);
  }
}

//####################################################################################################
// QSDcontract with index-based contraction scaling
//####################################################################################################

// Auxiliary function (1)
template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
double f_permuted_scale_axb
(const TVector<Q, NA>& a_permuted_qindex, const IVector<NA>& a_permute,
 const TVector<Q, NB>& b_permuted_qindex, const IVector<NB>& b_permute,
 const TVector<Q, NC>& c_qindex,
 const function<double(const TVector<Q, NA>&,
                       const TVector<Q, NB>&,
                       const TVector<Q, NC>&)>& f_scale)
{
  TVector<Q, NA> a_qindex;
  for(int i = 0; i < NA; ++i) a_qindex[a_permute[i]] = a_permuted_qindex[i];
  TVector<Q, NB> b_qindex;
  for(int i = 0; i < NB; ++i) b_qindex[b_permute[i]] = b_permuted_qindex[i];
  return f_scale(a_qindex, b_qindex, c_qindex);
}

// Auxiliary function (2)
template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
double f_permuted_scale_c
(const TVector<Q, NA>& a_qindex,
 const TVector<Q, NB>& b_qindex,
 const TVector<Q, NC>& c_permuted_qindex, const IVector<NC>& c_permute,
 const function<double(const TVector<Q, NA>&,
                       const TVector<Q, NB>&,
                       const TVector<Q, NC>&)>& f_scale)
{
  TVector<Q, NC> c_qindex;
  for(int i = 0; i < NC; ++i) c_qindex[c_permute[i]] = c_permuted_qindex[i];
  return f_scale(a_qindex, b_qindex, c_qindex);
}

template<size_t NA, size_t NB, size_t NC, size_t K, class Q = Quantum>
void QSDcontract
(const function<double(const TVector<Q, NA>&,
                       const TVector<Q, NB>&,
                       const TVector<Q, NC>&)>& f_scale,
 const double& alpha, const QSDArray<NA, Q>& a, const IVector<K>& a_contract,
                      const QSDArray<NB, Q>& b, const IVector<K>& b_contract,
 const double& beta,        QSDArray<NC, Q>& c)
{
  IVector<NA> a_permute;
  IVector<NB> b_permute;
  unsigned int jobs = get_contract_jobs(a.shape(), a_contract, a_permute,
                                        b.shape(), b_contract, b_permute);

  QSDArray<NA, Q> a_ref;
  if(jobs & JOBMASK_A_PMUTE) QSDpermute(a, a_permute, a_ref);
  else                       a_ref.reference(a);

  BTAS_TRANSPOSE TransA;
  if(jobs & JOBMASK_A_TRANS) TransA = Trans;
  else                       TransA = NoTrans;

  QSDArray<NB, Q> b_ref;
  if(jobs & JOBMASK_B_PMUTE) QSDpermute(b, b_permute, b_ref);
  else                       b_ref.reference(b);

  BTAS_TRANSPOSE TransB;
  if(jobs & JOBMASK_B_TRANS) TransB = Trans;
  else                       TransB = NoTrans;

  IVector<NC> c_permute;
  for(int i = 0; i < NC; ++i) c_permute[i] = i;
  function<double(const TVector<Q, NA>&, const TVector<Q, NB>&, const TVector<Q, NC>&)>
  q_scale = bind(f_permuted_scale_axb<NA, NB, NC>, _1, a_permute, _2, b_permute, _3, f_scale);

  switch(jobs & JOBMASK_BLAS_TYPE) {
    case(0):
      QSDgemv(q_scale, TransA, alpha, a_ref, b_ref, beta, c);
      break;
    case(1):
      QSDgemv(q_scale, TransB, alpha, b_ref, a_ref, beta, c);
      break;
    case(2):
      QSDgemm(q_scale, TransA, TransB, alpha, a_ref, b_ref, beta, c);
      break;
    default:
      BTAS_THROW(false, "btas::QSDcontract unknown blas job type returned");
  }
}

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void QSDindexed_contract
(const function<double(const TVector<Q, NA>&,
                       const TVector<Q, NB>&,
                       const TVector<Q, NC>&)>& f_scale,
 const double& alpha, const QSDArray<NA, Q>& a, const IVector<NA>& a_symbols,
                      const QSDArray<NB, Q>& b, const IVector<NB>& b_symbols,
 const double& beta,        QSDArray<NC, Q>& c, const IVector<NC>& c_symbols)
{
  const size_t K = (NA + NB - NC)/2;
  IVector<K> a_contract;
  IVector<K> b_contract;
  IVector<NC> axb_symbols;
  indexed_contract_shape(a_symbols, a_contract, b_symbols, b_contract, axb_symbols);

  if(c_symbols == axb_symbols) {
    QSDcontract(alpha, a, a_contract, b, b_contract, beta, c);
  }
  else {
    QSDArray<NC, Q> axb;
    IVector<NC> c_permute;
    indexed_permute_shape(c_symbols, axb_symbols, c_permute);
    if(c.size())
      QSDpermute(c, c_permute, axb);

    function<double(const TVector<Q, NA>&, const TVector<Q, NB>&, const TVector<Q, NC>&)>
    q_scale = bind(f_permuted_scale_c<NA, NB, NC>, _1, _2, _3, c_permute, f_scale);

    QSDcontract(q_scale, alpha, a, a_contract, b, b_contract, beta, axb);
    QSDindexed_permute(axb, axb_symbols, c, c_symbols);
  }
}

};

#endif // _BTAS_CXX11_QSDCONTRACT_H
