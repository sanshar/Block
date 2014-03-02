#ifndef _BTAS_CXX11_SDPERMUTE_H
#define _BTAS_CXX11_SDPERMUTE_H 1

#include <vector>
#include <set>
#include <algorithm>

#include <btas/TVector.h>
#include <btas/btas_permute_shape.h>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/Darglist.h>

#include <btas/SPARSE/SDArray.h>
#include <btas/SPARSE/SDblas.h>

namespace btas {

//! Call sparse-array permutation (serial run)
template<size_t N>
void serial_SDpermute
(const SDArray<N>& x, const IVector<N>& pindex, SDArray<N>& y)
{
  for(typename SDArray<N>::const_iterator ix = x.begin(); ix != x.end(); ++ix) {
    IVector<N> x_index = x.index(ix->first);
    IVector<N> y_index = permute(x_index, pindex);
    typename SDArray<N>::iterator iy = y.reserve(y_index);
    Dpermute(*(ix->second), pindex, *(iy->second));
  }
}

//! Call sparse-array permutation (threaded)
template<size_t N>
void thread_SDpermute
(const SDArray<N>& x, const IVector<N>& pindex, SDArray<N>& y)
{
  std::vector<DpermuteArglist<N>> task_list;
  task_list.reserve(x.nnz());
  for(typename SDArray<N>::const_iterator ix = x.begin(); ix != x.end(); ++ix) {
    IVector<N> x_index = x.index(ix->first);
    IVector<N> y_index = permute(x_index, pindex);
    typename SDArray<N>::iterator iy = y.reserve(y_index);
    task_list.push_back(DpermuteArglist<N>(ix->second, pindex, iy->second));
  }
  parallel_call(task_list);
}

//! Call sparse-array permutation (wrapper)
template<size_t N>
void SDpermute
(const SDArray<N>& x, const IVector<N>& pindex, SDArray<N>& y)
{
  std::set<int> iset(pindex.begin(), pindex.end());
  if(  iset.size()    != N ) BTAS_THROW(false, "btas::SDpermute: found duplicate index");
  if(*(iset.rbegin()) >= N ) BTAS_THROW(false, "btas::SDpermute: found out-of-range index");

  if(std::equal(iset.begin(), iset.end(), pindex.begin())) {
    SDcopy(x, y);
  }
  else {
    y.resize(permute(x.dshape(), pindex), false);
#ifdef _SERIAL
    serial_SDpermute(x, pindex, y);
#else
    thread_SDpermute(x, pindex, y);
#endif
  }
}

//! Call sparse-array indexed-permutation (wrapper)
template<size_t N>
void SDindexed_permute
(const SDArray<N>& x, const IVector<N>& x_symbols,
       SDArray<N>& y, const IVector<N>& y_symbols)
{
  if(x_symbols == y_symbols) {
    SDcopy(x, y);
  }
  else {
    IVector<N> pindex;
    indexed_permute_shape(x_symbols, y_symbols, pindex);
    SDpermute(x, pindex, y);
  }
}

}; // namespace btas

#endif // _BTAS_CXX11_SDPERMUTE_H
