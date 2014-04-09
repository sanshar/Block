#ifndef _BTAS_CXX11_DDIAGONAL_H
#define _BTAS_CXX11_DDIAGONAL_H 1

#include <set>

#include <btas/btas.h>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/Dreindex.h>

namespace btas {

//! Function to get diagonal elements
/*! Note that the order of d_index affects to the index of 'b'
 *  a(i,j,k,l) with d_index(j,l) returns b(i,j,k) = a(i,j,k,j)
 *  a(i,j,k,l) with d_index(l,j) returns b(i,k,l) = a(i,l,k,l)
 */
template<size_t NA, size_t K>
void Ddiagonal
(const DArray<NA>& a, const IVector<K>& d_index, DArray<NA-K+1>& b) {
  const size_t NB = NA - K + 1;
  // check consistency b/w 'a' shape and diag. index
  if(a.size() == 0)
    BTAS_THROW(false, "btas::Ddiagonal: array data not found");
  const IVector<NA>& a_shape  = a.shape ();
  const IVector<NA>& a_stride = a.stride();
  for(int i = 1; i < K; ++i) {
    if(a_shape[d_index[0]] != a_shape[d_index[i]])
      BTAS_THROW(false, "btas::Ddiagonal: diagonal indices must have the same size");
  }
  // calc. unique index
  IVector<NB> u_index;
  std::set<int> iset(d_index.begin(), d_index.end());
  int id = d_index[0];
  int nb = 0;
  // search unique index before id
  for(int i = 0; i < id; ++i)
    if(iset.find(i) == iset.end()) u_index[nb++] = i;
  // diagonal index
  u_index[nb++] = id;
  // search unique index after id
  for(int i = id+1; i < NA; ++i)
    if(iset.find(i) == iset.end()) u_index[nb++] = i;
  // 'b' shape and stride to pick diag. elements from 'a'
  IVector<NB> b_shape;
  IVector<NB> diagstr;
  for(int i = 0; i < NB; ++i) {
    b_shape[i] = a_shape [u_index[i]];
    diagstr[i] = a_stride[u_index[i]];
    if(u_index[i] == id)
      for(int j = 1; j < K; ++j) diagstr[i] += a_stride[d_index[j]];
  }
  // call Dreindex by diag. stride to pick diag. elements up
  b.resize(b_shape);
  Dreindex(a.data(), b.data(), diagstr, b_shape);
}

}; // namespace btas

#endif // _BTAS_CXX11_DIAGONAL_H
