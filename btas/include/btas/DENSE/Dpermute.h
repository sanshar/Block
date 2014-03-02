#ifndef _BTAS_CXX11_DPERMUTE_H
#define _BTAS_CXX11_DPERMUTE_H 1

#include <set>

#include <btas/btas.h>
#include <btas/btas_permute_shape.h>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/Dreindex.h>

namespace btas {

//! Permutation double precision real dense array
template<size_t N>
void Dpermute
(const DArray<N>& x, const IVector<N>& pindx, DArray<N>& y)
{
  if(x.size() == 0)
    BTAS_THROW(false, "btas::Dpermute: array data not found");

  std::set<int> iset(pindx.begin(), pindx.end());

  if(  iset.size()    != N )
    BTAS_THROW(false, "btas::Dpermute: duplicate indices found");
  if(*(iset.rbegin()) >= N )
    BTAS_THROW(false, "btas::Dpermute: out-of-range indices found");

  if(std::equal(iset.begin(), iset.end(), pindx.begin())) {
    Dcopy(x, y);
  }
  else {
    IVector<N> xstr;
    IVector<N> yshape;
    permute_shape(pindx, x.shape(), xstr, yshape);
    y.resize(yshape);
    Dreindex(x.data(), y.data(), xstr, yshape);
  }
}

//! Indexed permutation for double precision real dense array
template<size_t N>
void Dindexed_permute
(const DArray<N>& x, const IVector<N>& x_symbols, DArray<N>& y, const IVector<N>& y_symbols)
{
  if(x_symbols == y_symbols) {
    Dcopy(x, y);
  }
  else {
    IVector<N> pindx;
    indexed_permute_shape(x_symbols, y_symbols, pindx);
    Dpermute(x, pindx, y);
  }
}

}; // namespace btas

#endif // _BTAS_CXX11_DPERMUTE_H
