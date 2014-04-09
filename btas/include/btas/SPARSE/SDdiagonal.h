#ifndef _BTAS_CXX11_SDDIAGONAL_H
#define _BTAS_CXX11_SDDIAGONAL_H 1

#include <set>

#include <btas/btas.h>
#include <btas/TVector.h>

#include <btas/DENSE/Ddiagonal.h>

#include <btas/SPARSE/SDArray.h>

namespace btas {

//! Function to get diagonal elements for sparse array
/*! Note that the order of d_index affects to the index of 'b'
 *  a(i,j,k,l) with d_index(j,l) returns b(i,j,k) = a(i,j,k,j)
 *  a(i,j,k,l) with d_index(l,j) returns b(i,k,l) = a(i,l,k,l)
 */
template<size_t NA, size_t K>
void SDdiagonal
(const SDArray<NA>& a, const IVector<K>& d_index, SDArray<NA-K+1>& b) {
  const size_t NB = NA - K + 1;
  if(a.size() == 0)
    BTAS_THROW(false, "btas::SDdiagonal: array data not found");
  // check consistency b/w 'a' shape and diag. index 'd_index'
  const TVector<Dshapes, NA>& a_dn_shape = a.dshape();
  for(int i = 1; i < K; ++i) {
    if(a_dn_shape[d_index[0]] != a_dn_shape[d_index[i]])
      BTAS_THROW(false, "btas::SDdiagonal: diagonal indices must have the same size");
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
  // calc. 'b' shape and resize
  TVector<Dshapes, NB> b_dn_shape;
  for(int i = 0; i < NB; ++i) b_dn_shape[i] = a_dn_shape[u_index[i]];
  b.resize(b_dn_shape, false);
  // pick diagonal elements
  for(typename SDArray<NA>::const_iterator ia = a.begin(); ia != a.end(); ++ia) {
    IVector<NA> a_index = a.index(ia->first);
    // check whether a_index is diagonal block
    bool is_diag_block = true;
    for(int i = 1; i < K; ++i) {
      if(a_index[id] != a_index[d_index[i]]) {
        is_diag_block = false; break;
      }
    }
    if(!is_diag_block) continue;
    // diag. block
    IVector<NB> b_index;
    for(int i = 0; i < NB; ++i) b_index[i] = a_index[u_index[i]];
    typename SDArray<NB>::iterator ib = b.reserve(b_index);
    if(ib == b.end())
      BTAS_THROW(false, "btas::SDdiagonal: requested block must be zero, could not be reserved");
    // call for dense array
    Ddiagonal((*ia->second), d_index, (*ib->second));
  }
}

}; // namespace btas

#endif // _BTAS_CXX11_SDIAGONAL_H
