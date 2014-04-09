#ifndef _BTAS_CXX11_DRIVER_PERMUTE_H
#define _BTAS_CXX11_DRIVER_PERMUTE_H 1

#include <btas/DENSE/Dpermute.h>
#include <btas/SPARSE/SDpermute.h>
#include <btas/QSPARSE/QSDpermute.h>

namespace btas {

template<size_t N>
inline void permute
(const DArray<N>& x, const IVector<N>& index, DArray<N>& y)
{
  Dpermute(x, index, y);
}

template<size_t N>
inline void indexed_permute
(const DArray<N>& x, const IVector<N>& x_symbols,
       DArray<N>& y, const IVector<N>& y_symbols)
{
  Dindexed_permute(x, x_symbols, y, y_symbols);
}

template<size_t N>
inline void permute
(const SDArray<N>& x, const IVector<N>& index, SDArray<N>& y)
{
  SDpermute(x, index, y);
}

template<size_t N>
inline void indexed_permute
(const SDArray<N>& x, const IVector<N>& x_symbols,
       SDArray<N>& y, const IVector<N>& y_symbols)
{
  SDindexed_permute(x, x_symbols, y, y_symbols);
}

template<size_t N, class Q = Quantum>
inline void permute
(const QSDArray<N, Q>& x, const IVector<N>& index, QSDArray<N, Q>& y)
{
  QSDpermute(x, index, y);
}

template<size_t N, class Q = Quantum>
inline void indexed_permute
(const QSDArray<N, Q>& x, const IVector<N>& x_symbols,
       QSDArray<N, Q>& y, const IVector<N>& y_symbols)
{
  QSDindexed_permute(x, x_symbols, y, y_symbols);
}

};

#endif // _BTAS_CXX11_DRIVER_PERMUTE_H
