#ifndef _BTAS_CXX11_DRIVER_CONTRACT_H
#define _BTAS_CXX11_DRIVER_CONTRACT_H 1

#include <btas/DENSE/Dcontract.h>
#include <btas/SPARSE/SDcontract.h>
#include <btas/QSPARSE/QSDcontract.h>

namespace btas {

template<size_t NA, size_t NB, size_t K>
inline void contract
(const double& alpha, const DArray<NA>& a, const IVector<K>& a_contract,
                      const DArray<NB>& b, const IVector<K>& b_contract,
 const double& beta,        DArray<NA+NB-K-K>& c)
{
  Dcontract(alpha, a, a_contract, b, b_contract, beta, c);
}

template<size_t NA, size_t NB, size_t NC>
inline void indexed_contract
(const double& alpha, const DArray<NA>& a, const IVector<NA>& a_symbols,
                      const DArray<NB>& b, const IVector<NB>& b_symbols,
 const double& beta,        DArray<NC>& c, const IVector<NC>& c_symbols)
{
  Dindexed_contract(alpha, a, a_symbols, b, b_symbols, beta, c, c_symbols);
}

template<size_t NA, size_t NB, size_t K>
inline void contract
(const double& alpha, const SDArray<NA>& a, const IVector<K>& a_contract,
                      const SDArray<NB>& b, const IVector<K>& b_contract,
 const double& beta,        SDArray<NA+NB-K-K>& c)
{
  SDcontract(alpha, a, a_contract, b, b_contract, beta, c);
}

template<size_t NA, size_t NB, size_t NC>
inline void indexed_contract
(const double& alpha, const SDArray<NA>& a, const IVector<NA>& a_symbols,
                      const SDArray<NB>& b, const IVector<NB>& b_symbols,
 const double& beta,        SDArray<NC>& c, const IVector<NC>& c_symbols)
{
  SDindexed_contract(alpha, a, a_symbols, b, b_symbols, beta, c, c_symbols);
}

template<size_t NA, size_t NB, size_t K, class Q = Quantum>
inline void contract
(const double& alpha, const QSDArray<NA, Q>& a, const IVector<K>& a_contract,
                      const QSDArray<NB, Q>& b, const IVector<K>& b_contract,
 const double& beta,        QSDArray<NA+NB-K-K, Q>& c)
{
  QSDcontract(alpha, a, a_contract, b, b_contract, beta, c);
}

template<size_t NA, size_t NB, size_t NC, class Q = Quantum>
inline void indexed_contract
(const double& alpha, const QSDArray<NA, Q>& a, const IVector<NA>& a_symbols,
                      const QSDArray<NB, Q>& b, const IVector<NB>& b_symbols,
 const double& beta,        QSDArray<NC, Q>& c, const IVector<NC>& c_symbols)
{
  QSDindexed_contract(alpha, a, a_symbols, b, b_symbols, beta, c, c_symbols);
}

};

#endif // _BTAS_CXX11_DRIVER_CONTRACT_H
