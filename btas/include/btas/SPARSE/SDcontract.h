#ifndef _BTAS_CXX11_SDCONTRACT_H
#define _BTAS_CXX11_SDCONTRACT_H 1

#include <btas/btas_contract_shape.h>

#include <btas/SPARSE/SDArray.h>
#include <btas/SPARSE/SDblas.h>
#include <btas/SPARSE/SDpermute.h>

namespace btas
{

//! Convenient sparse-array contraction function
template<size_t NA, size_t NB, size_t K>
void SDcontract
(const double& alpha,
 const SDArray<NA>& a, const IVector<K>& a_contract,
 const SDArray<NB>& b, const IVector<K>& b_contract,
 const double& beta,
       SDArray<NA+NB-K-K>& c)
{
  IVector<NA> a_pindex;
  IVector<NB> b_pindex;
  unsigned int jobs = get_contract_jobs(a.shape(), a_contract, a_pindex, b.shape(), b_contract, b_pindex);

  SDArray<NA> a_ref;
  if(jobs & JOBMASK_A_PMUTE) SDpermute(a, a_pindex, a_ref);
  else                       a_ref.reference(a);

  BTAS_TRANSPOSE TransA;
  if(jobs & JOBMASK_A_TRANS) TransA = Trans;
  else                       TransA = NoTrans;

  SDArray<NB> b_ref;
  if(jobs & JOBMASK_B_PMUTE) SDpermute(b, b_pindex, b_ref);
  else                       b_ref.reference(b);

  BTAS_TRANSPOSE TransB;
  if(jobs & JOBMASK_B_TRANS) TransB = Trans;
  else                       TransB = NoTrans;

  switch(jobs & JOBMASK_BLAS_TYPE) {
    case(0):
      SDgemv(TransA, alpha, a_ref, b_ref, beta, c);
      break;
    case(1):
      SDgemv(TransB, alpha, b_ref, a_ref, beta, c);
      break;
    case(2):
      SDgemm(TransA, TransB, alpha, a_ref, b_ref, beta, c);
      break;
    default:
      BTAS_THROW(false, "btas::SDcontract unknown blas job type returned");
  }
}

//! Convenient sparse-array indexed-contraction function
template<size_t NA, size_t NB, size_t NC>
void SDindexed_contract
(const double& alpha,
 const SDArray<NA>& a, const IVector<NA>& a_symbols,
 const SDArray<NB>& b, const IVector<NB>& b_symbols,
 const double& beta,
       SDArray<NC>& c, const IVector<NC>& c_symbols)
{
  const size_t K = (NA + NB - NC)/2;
  IVector<K> a_contract;
  IVector<K> b_contract;
  IVector<NC> axb_symbols;
  indexed_contract_shape(a_symbols, a_contract, b_symbols, b_contract, axb_symbols);

  if(c_symbols == axb_symbols) {
    SDcontract(alpha, a, a_contract, b, b_contract, beta, c);
  }
  else {
    SDArray<NC> axb;
    // permute if necessary
    if(c.size() > 0)
      SDindexed_permute(c, c_symbols, axb, axb_symbols);
    // array contraction
    SDcontract(alpha, a, a_contract, b, b_contract, beta, axb);
    // permute back to c_symbols order
    SDindexed_permute(axb, axb_symbols, c, c_symbols);
  }
}

}; // namespace btas

#endif // _BTAS_CXX11_SDCONTRACT_H
