
//####################################################################################################
/*! \file  btas_contract_qshape.h
 *  \brief Quantum number contraction for quantum number-based sparse BLAS calls */
//####################################################################################################

#ifndef _BTAS_CXX11_CONTRACT_QSHAPE_H
#define _BTAS_CXX11_CONTRACT_QSHAPE_H 1

#include <btas/btas.h>
#include <btas/btas_contract_shape.h>

#include <btas/QSPARSE/Qshapes.h>

namespace btas {

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void gemv_contract_qshape
(const BTAS_TRANSPOSE& TransA,
 const Q& a_qnum, const TVector<Qshapes<Q>, NA>& a_qshape,
 const Q& b_qnum, const TVector<Qshapes<Q>, NB>& b_qshape,
       Q& c_qnum,       TVector<Qshapes<Q>, NC>& c_qshape)
{
  c_qnum = Q::zero();
  if(TransA == NoTrans) {
    c_qnum =  a_qnum * b_qnum;
    for(int i = 0; i < NC; ++i) c_qshape[i] =  a_qshape[i];
    for(int i = 0; i < NB; ++i)
      if(a_qshape[i+NC] != -b_qshape[i])
        BTAS_THROW(false, "btas::gemv_contract_qshape contraction of quantum numbers failed");
  }
  else if(TransA == ConjTrans) {
    c_qnum = -a_qnum * b_qnum;
    for(int i = 0; i < NC; ++i) c_qshape[i] = -a_qshape[i+NB];
    for(int i = 0; i < NB; ++i)
      if(a_qshape[i]    !=  b_qshape[i])
        BTAS_THROW(false, "btas::gemv_contract_qshape contraction of quantum numbers failed");
  }
  else {
    c_qnum =  a_qnum * b_qnum;
    for(int i = 0; i < NC; ++i) c_qshape[i] =  a_qshape[i+NB];
    for(int i = 0; i < NB; ++i)
      if(a_qshape[i]    != -b_qshape[i])
        BTAS_THROW(false, "btas::gemv_contract_qshape contraction of quantum numbers failed");
  }
}

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void ger_contract_qshape
(const Q& a_qnum, const TVector<Qshapes<Q>, NA>& a_qshape,
 const Q& b_qnum, const TVector<Qshapes<Q>, NB>& b_qshape,
       Q& c_qnum,       TVector<Qshapes<Q>, NC>& c_qshape)
{
  c_qnum = a_qnum * b_qnum;
  for(int i = 0; i < NA; ++i) c_qshape[i]    = a_qshape[i];
  for(int i = 0; i < NB; ++i) c_qshape[i+NA] = b_qshape[i];
}

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
void gemm_contract_qshape
(const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
 const Q& a_qnum, const TVector<Qshapes<Q>, NA>& a_qshape,
 const Q& b_qnum, const TVector<Qshapes<Q>, NB>& b_qshape,
       Q& c_qnum,       TVector<Qshapes<Q>, NC>& c_qshape)
{
  const size_t K = (NA + NB - NC)/2;

  TVector<Qshapes<Q>, K> k_qshape;
  c_qnum = Q::zero();

  if(TransA == NoTrans) {
    c_qnum =  a_qnum;
    for(int i = 0; i < NA - K; ++i) c_qshape[i] =  a_qshape[i];
    for(int i = 0; i < K;      ++i) k_qshape[i] =  a_qshape[i+NA-K];
  }
  else if(TransA == ConjTrans) {
    c_qnum = -a_qnum;
    for(int i = 0; i < NA - K; ++i) c_qshape[i] = -a_qshape[i+K];
    for(int i = 0; i < K;      ++i) k_qshape[i] = -a_qshape[i];
  }
  else {
    c_qnum =  a_qnum;
    for(int i = 0; i < NA - K; ++i) c_qshape[i] =  a_qshape[i+K];
    for(int i = 0; i < K;      ++i) k_qshape[i] =  a_qshape[i];
  }

  if(TransB == NoTrans) {
    c_qnum =  b_qnum * c_qnum;
    for(int i = 0; i < NB - K; ++i) c_qshape[i+NA-K] =  b_qshape[i+K];
    for(int i = 0; i < K; ++i)
      if(k_qshape[i] != -b_qshape[i])
        BTAS_THROW(false, "btas::gemm_contract_qshape contraction of quantum numbers failed");
  }
  else if(TransB == ConjTrans) {
    c_qnum = -b_qnum * c_qnum;
    for(int i = 0; i < NB - K; ++i) c_qshape[i+NA-K] = -b_qshape[i];
    for(int i = 0; i < K; ++i)
      if(k_qshape[i] !=  b_qshape[i+NB-K])
        BTAS_THROW(false, "btas::gemm_contract_qshape contraction of quantum numbers failed");
  }
  else {
    c_qnum =  b_qnum * c_qnum;
    for(int i = 0; i < NB - K; ++i) c_qshape[i+NA-K] =  b_qshape[i];
    for(int i = 0; i < K; ++i)
      if(k_qshape[i] != -b_qshape[i+NB-K])
        BTAS_THROW(false, "btas::gemm_contract_qshape contraction of quantum numbers failed");
  }
}

//####################################################################################################
// Index-based contraction scaling function
//####################################################################################################

//! Wrapper function to return Clebsch-Gordan coefficient depends on block index
template<typename T, size_t NA, size_t NB, size_t NC, class Q = Quantum>
double f_indexbase_scale
(const function<T(const TVector<Q, NA>&,
                  const TVector<Q, NB>&,
                  const TVector<Q, NC>&)>& f_scale,
 const TVector<Qshapes<Q>, NA>& a_qshape, const IVector<NA>& a_index,
 const TVector<Qshapes<Q>, NB>& b_qshape, const IVector<NB>& b_index,
 const TVector<Qshapes<Q>, NC>& c_qshape, const IVector<NC>& c_index)
{
  TVector<Q, NA> a_qindex;
  for(int i = 0; i < NA; ++i) a_qindex[i] = a_qshape[i][a_index[i]];
  TVector<Q, NB> b_qindex;
  for(int i = 0; i < NB; ++i) b_qindex[i] = b_qshape[i][b_index[i]];
  TVector<Q, NC> c_qindex;
  for(int i = 0; i < NC; ++i) c_qindex[i] = c_qshape[i][c_index[i]];
  return f_scale(a_qindex, b_qindex, c_qindex);
}

}; // namespace btas

#endif // _BTAS_CXX11_CONTRACT_QSHAPE_H
