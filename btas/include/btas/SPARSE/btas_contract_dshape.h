
//####################################################################################################
/*! \file  btas_contract_dshape.h
 *  \brief Block-shape contraction for sparse array BLAS calls */
//####################################################################################################

#ifndef _BTAS_CXX11_CONTRACT_DSHAPE_H
#define _BTAS_CXX11_CONTRACT_DSHAPE_H 1

#include <btas/btas.h>
#include <btas/btas_contract_shape.h>

namespace btas {

//! Compare Dshapes ignoring 0-sized index
inline bool _btas_contract_allowed(const Dshapes& x, const Dshapes& y) {
  bool _equal = (x.size() == y.size());
  for(size_t i = 0; _equal && i < x.size(); ++i)
    if(x[i] > 0 && y[i] > 0) _equal &= (x[i] == y[i]);
  return _equal;
}

template<size_t NA, size_t NB, size_t NC>
void gemv_contract_dshape
(const BTAS_TRANSPOSE& TransA,
 const TVector<Dshapes, NA>& a_dshape,
 const TVector<Dshapes, NB>& b_dshape,
       TVector<Dshapes, NC>& c_dshape)
{
  if(TransA == NoTrans) {
    for(int i = 0; i < NC; ++i) c_dshape[i] = a_dshape[i];
    for(int i = 0; i < NB; ++i)
      BTAS_THROW(_btas_contract_allowed(a_dshape[i+NC], b_dshape[i]), "btas::gemv_contract_dshape contraction of dense-array shape failed");
  }
  else {
    for(int i = 0; i < NC; ++i) c_dshape[i] = a_dshape[i+NB];
    for(int i = 0; i < NB; ++i)
      BTAS_THROW(_btas_contract_allowed(a_dshape[i],    b_dshape[i]), "btas::gemv_contract_dshape contraction of dense-array shape failed");
  }
}

template<size_t NA, size_t NB, size_t NC>
void ger_contract_dshape
(const TVector<Dshapes, NA>& a_dshape,
 const TVector<Dshapes, NB>& b_dshape,
       TVector<Dshapes, NC>& c_dshape)
{
  for(int i = 0; i < NA; ++i) c_dshape[i]    = a_dshape[i];
  for(int i = 0; i < NB; ++i) c_dshape[i+NA] = b_dshape[i];
}

template<size_t NA, size_t NB, size_t NC>
void gemm_contract_dshape
(const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
 const TVector<Dshapes, NA>& a_dshape,
 const TVector<Dshapes, NB>& b_dshape,
       TVector<Dshapes, NC>& c_dshape)
{
  const size_t K = (NA + NB - NC)/2;

  TVector<Dshapes, K> k_dshape;

  if(TransA == NoTrans) {
    for(int i = 0; i < NA-K; ++i) c_dshape[i] = a_dshape[i];
    for(int i = 0; i < K;    ++i) k_dshape[i] = a_dshape[i+NA-K];
  }
  else {
    for(int i = 0; i < NA-K; ++i) c_dshape[i] = a_dshape[i+K];
    for(int i = 0; i < K;    ++i) k_dshape[i] = a_dshape[i];
  }

  if(TransB == NoTrans) {
    for(int i = 0; i < NB-K; ++i) c_dshape[i+NA-K] = b_dshape[i+K];
    for(int i = 0; i < K; ++i)
      BTAS_THROW(_btas_contract_allowed(k_dshape[i], b_dshape[i]),      "btas::gemm_contract_dshape contraction of dense-array shape failed");
  }
  else {
    for(int i = 0; i < NB-K; ++i) c_dshape[i+NA-K] = b_dshape[i];
    for(int i = 0; i < K; ++i)
      BTAS_THROW(_btas_contract_allowed(k_dshape[i], b_dshape[i+NB-K]), "btas::gemm_contract_dshape contraction of dense-array shape failed");
  }
}

}; // namespace btas

#endif // _BTAS_CXX11_CONTRACT_DSHAPE_H
