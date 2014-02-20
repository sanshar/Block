
//
/*! \file  QSDpermute.h
 *  \brief Array permutation for QSDArray
 */

#ifndef _BTAS_CXX11_QSDPERMUTE_H
#define _BTAS_CXX11_QSDPERMUTE_H 1

#include <set>
#include <algorithm>

#include <btas/SPARSE/SDpermute.h>

#include <btas/QSPARSE/QSDArray.h>

namespace btas {

template<size_t N, class Q = Quantum>
void QSDpermute
(const QSDArray<N, Q>& x, const IVector<N>& pindex, QSDArray<N, Q>& y)
{
  std::set<int> iset(pindex.begin(), pindex.end());
  if(  iset.size()    != N ) BTAS_THROW(false, "btas::QSDpermute: found duplicate index");
  if(*(iset.rbegin()) >= N ) BTAS_THROW(false, "btas::QSDpermute: found out-of-range index");

  if(std::equal(iset.begin(), iset.end(), pindex.begin())) {
    QSDcopy(x, y);
  }
  else {
    y.resize(x.q(), permute(x.qshape(), pindex), permute(x.dshape(), pindex), false);
#ifdef _SERIAL
    serial_SDpermute(x, pindex, y);
#else
    thread_SDpermute(x, pindex, y);
#endif
  }
}

template<size_t N, class Q = Quantum>
void QSDindexed_permute
(const QSDArray<N, Q>& x, const IVector<N>& x_symbols,
       QSDArray<N, Q>& y, const IVector<N>& y_symbols)
{
  if(x_symbols == y_symbols) {
    QSDcopy(x, y);
  }
  else {
    IVector<N> pindex;
    indexed_permute_shape(x_symbols, y_symbols, pindex);
    QSDpermute(x, pindex, y);
  }
}

}; // namespace btas

#endif // _BTAS_CXX11_QSDPERMUTE_H
