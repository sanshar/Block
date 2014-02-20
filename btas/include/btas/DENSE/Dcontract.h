#ifndef _BTAS_CXX11_DCONTRACT_H
#define _BTAS_CXX11_DCONTRACT_H 1

#include <btas/btas.h>
#include <btas/btas_contract_shape.h>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/Dblas.h>
#include <btas/DENSE/Dpermute.h>

namespace btas {

//! Convenient contraction function for DArray
template<size_t NA, size_t NB, size_t K>
void Dcontract
(const double& alpha,
 const DArray<NA>& a, const IVector<K>& a_contract,
 const DArray<NB>& b, const IVector<K>& b_contract,
 const double& beta,
       DArray<NA+NB-K-K>& c) {

  IVector<NA> a_pindex;
  IVector<NB> b_pindex;
  unsigned int jobs = get_contract_jobs(a.shape(), a_contract, a_pindex, b.shape(), b_contract, b_pindex);

  DArray<NA> a_ref;
  if(jobs & JOBMASK_A_PMUTE) Dpermute(a, a_pindex, a_ref);
  else                       a_ref.reference(a);

  BTAS_TRANSPOSE transa;
  if(jobs & JOBMASK_A_TRANS) transa = Trans;
  else                       transa = NoTrans;

  DArray<NB> b_ref;
  if(jobs & JOBMASK_B_PMUTE) Dpermute(b, b_pindex, b_ref);
  else                       b_ref.reference(b);

  BTAS_TRANSPOSE transb;
  if(jobs & JOBMASK_B_TRANS) transb = Trans;
  else                       transb = NoTrans;

  switch(jobs & JOBMASK_BLAS_TYPE) {
    case(0):
      Dgemv(transa, alpha, a_ref, b_ref, beta, c);
      break;
    case(1):
      Dgemv(transb, alpha, b_ref, a_ref, beta, c);
      break;
    case(2):
      Dgemm(transa, transb, alpha, a_ref, b_ref, beta, c);
      break;
    default:
      BTAS_THROW(false, "btas::Dcontract: unknown BLAS job type returned");
  }
}

//! Convenient indexed contraction function for DArray
template<size_t NA, size_t NB, size_t NC>
void Dindexed_contract
(const double& alpha,
 const DArray<NA>& a, const IVector<NA>& a_symbols,
 const DArray<NB>& b, const IVector<NB>& b_symbols,
 const double& beta,
       DArray<NC>& c, const IVector<NC>& c_symbols) {

  const size_t K = (NA + NB - NC)/2;

  IVector<K> a_contract;
  IVector<K> b_contract;
  IVector<NC> axb_symbols;
  indexed_contract_shape(a_symbols, a_contract, b_symbols, b_contract, axb_symbols);

  if(c_symbols == axb_symbols) {
    Dcontract(alpha, a, a_contract, b, b_contract, beta, c);
  }
  else {
    DArray<NC> axb;
    if(c.size() > 0) {
      Dindexed_permute(c, c_symbols, axb, axb_symbols);
    }
    Dcontract(alpha, a, a_contract, b, b_contract, beta, axb);
    Dindexed_permute(axb, axb_symbols, c, c_symbols);
  }
}

}; // namespace btas

#endif // _BTAS_CXX11_DCONTRACT_H
