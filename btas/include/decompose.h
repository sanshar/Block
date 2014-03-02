
//
/*! \file  decompose.h
 *  \brief BTAS driver routines for tensor decomposition.
 *
 *  - Provide convenient wrappers to tensor factorization using SVD.
 *  - Arrow direction in QSDgesvd is always taken to be LeftArray
 *  - Truncated singular values are normalized, thus the norm of factorized tensor is preserved.
 */

#ifndef _BTAS_CXX11_DRIVER_DECOMPOSE_H
#define _BTAS_CXX11_DRIVER_DECOMPOSE_H 1

#include <set>
#include <map>

#include <btas/DENSE/Dlapack.h>
#include <btas/DENSE/Dpermute.h>

#include <btas/QSPARSE/QSDlapack.h>
#include <btas/QSPARSE/QSDpermute.h>

namespace btas {

template<size_t N, size_t K>
void decompose_shape
(const IVector<K>& index, IVector<N>& a_permute)
{
  std::set<int> iset(index.begin(), index.end());
  if(iset.size() < K)
    BTAS_THROW(false, "btas::decompose_shape: duplicate indices found in index");
  if((*iset.begin()) < 0 || (*iset.rbegin()) >= N)
    BTAS_THROW(false, "btas::decompose_shape: index is out-of-range");
  int n = 0;
  for(int i = 0; i < N; ++i)
    if(iset.find(i) == iset.end())
      a_permute[n++] = i;
  for(int i = 0; i < K; ++i)
      a_permute[i+N-K] = index[i];
}

template<size_t NA, size_t NL>
void indexed_decompose_shape
(const IVector<NA>& a_symbols, const IVector<NL>& l_symbols, const IVector<NA-NL+2>& r_symbols,
       IVector<NA>& a_permute,       IVector<NL>& l_permute,       IVector<NA-NL+2>& r_permute)
{
  const size_t NR = NA - NL + 2;
  const size_t K  = NR - 1;

  std::set <int> aset(a_symbols.begin(), a_symbols.end());

  std::pair<int, int> lshr;
  std::map <int, int> lmap;
  int nl = 0;
  for(int i = 0; i < NL; ++i) {
    if(aset.find(l_symbols[i]) != aset.end())
      lmap.insert(std::make_pair(l_symbols[i], nl++));
    else
      lshr = std::make_pair(l_symbols[i], i);
  }
  if(nl != lmap.size())
    BTAS_THROW(false, "btas::indexed_decompose_shape: duplicate indices found in l_symbols");
  if(nl != NL-1)
    BTAS_THROW(false, "btas::indexed_decompose_shape: there's missing index in l_symbols");

  std::pair<int, int> rshr;
  std::map <int, int> rmap;
  int nr = 0;
  for(int i = 0; i < NR; ++i) {
    if(aset.find(r_symbols[i]) != aset.end()) {
      if(lmap.find(r_symbols[i]) != lmap.end())
        BTAS_THROW(false, "btas::indexed_decompose_shape: found (new) factorized index in a_symbols");
      rmap.insert(std::make_pair(r_symbols[i], nr++));
    }
    else {
      rshr = std::make_pair(r_symbols[i], i);
    }
  }
  if(nr != rmap.size())
    BTAS_THROW(false, "btas::indexed_decompose_shape: duplicate indices found in r_symbols");
  if(nr != NR-1)
    BTAS_THROW(false, "btas::indexed_decompose_shape: there's missing index in r_symbols");

  if(lshr.first != rshr.first)
    BTAS_THROW(false, "btas::indexed_decompose_shape: there's no shared index between l_symbols and r_symbols");

  int n = 0;
  for(int i = 0; i < NA; ++i) {
    typename std::map<int, int>::iterator ilm = lmap.find(a_symbols[i]);
    if(ilm != lmap.end()) a_permute[ilm->second]      = i;
    typename std::map<int, int>::iterator irm = rmap.find(a_symbols[i]);
    if(irm != rmap.end()) a_permute[irm->second+NL-1] = i;
  }

  nl = lshr.second;
  l_permute[nl] = NL-1;
  for(int i = 0;      i < nl;   ++i) l_permute[i] = i;
  for(int i = nl + 1; i < NL-1; ++i) l_permute[i] = i-1;

  nr = rshr.second;
  r_permute[nr] = 0;
  for(int i = 0;      i < nr;   ++i) r_permute[i] = i+1;
  for(int i = nr + 1; i < NR-1; ++i) r_permute[i] = i;
}

// This should not be called directly
template<size_t N>
void _lsv_truncate(const DArray<N>& lsv, DArray<N>& lsv_tr, int D = 0)
{
  if(D == 0) return;
  IVector<N> ltr_shape(lsv.shape());
  ltr_shape[N-1] = D;
  lsv_tr.resize(ltr_shape);
  IVector<N> l_lbound = make_vector<N>(0);
  IVector<N> l_ubound;
  for(int i = 0; i < N; ++i) l_ubound[i] = ltr_shape[i]-1;
  copy(DSubArray<N>(lsv, l_lbound, l_ubound), lsv_tr);
}

// This should not be called directly
template<size_t N>
void _rsv_truncate(const DArray<N>& rsv, DArray<N>& rsv_tr, int D = 0)
{
  if(D == 0) return;
  IVector<N> rtr_shape(rsv.shape());
  rtr_shape[0] = D;
  rsv_tr.resize(rtr_shape);
  IVector<N> l_lbound = make_vector<N>(0);
  IVector<N> l_ubound;
  for(int i = 0; i < N; ++i) l_ubound[i] = rtr_shape[i]-1;
  copy(DSubArray<N>(rsv, l_lbound, l_ubound), rsv_tr);
}

template<size_t NA, size_t K>
double decompose
(const DArray<NA>& a, const IVector<K>& index, DArray<NA-K+1>& l, DArray<K+1>& r, int D = 0)
{
  const size_t NL = NA - K + 1;
  const size_t NR = K + 1;

  IVector<NA> a_permute;
  decompose_shape(index, a_permute);

  DArray<NA> a_copy;
  Dpermute(a, a_permute, a_copy);

  DArray<1>  s;
  DArray<NL> l_tmp;
  DArray<NL> r_tmp;
  Dgesvd(a_copy, s, l_tmp, r_tmp);
  Ddidm(s, r_tmp);
  _lsv_truncate(l_tmp, l, D);
  _rsv_truncate(r_tmp, r, D);
  // Selected norm
  double snorm = 0.0;
  for(int i = 0; i < D; ++i)        snorm += s(i) * s(i);
  // Discarded norm
  double dnorm = 0.0;
  for(int i = D; i < s.size(); ++i) dnorm += s(i) * s(i);
  // Return weight
  return (snorm / (snorm + dnorm));
}

template<size_t NA, size_t NL>
double indexed_decompose
(const DArray<NA>&      a, const IVector<NA>&      a_symbols,
       DArray<NL>&      l, const IVector<NL>&      l_symbols,
       DArray<NA-NL+2>& r, const IVector<NA-NL+2>& r_symbols, int D = 0)
{
  const size_t NR = NA - NL + 2;

  IVector<NA> a_permute;
  IVector<NL> l_permute;
  IVector<NR> r_permute;
  indexed_decompose_shape(a_symbols, l_symbols, r_symbols,
                          a_permute, l_permute, r_permute);

  DArray<NA> a_copy;
  Dpermute(a, a_permute, a_copy);

  DArray<1>  s;
  Dgesvd(a_copy, s, l, r);
  Ddidm(s, r);
  DArray<NL> l_tmp;
  DArray<NL> r_tmp;
  _lsv_truncate(l, l_tmp, D);
  _rsv_truncate(r, r_tmp, D);

  Dpermute(l_tmp, l_permute, l);
  Dpermute(r_tmp, r_permute, r);

  // Selected norm
  double snorm = 0.0;
  for(int i = 0; i < D; ++i)        snorm += s(i) * s(i);
  // Discarded norm
  double dnorm = 0.0;
  for(int i = D; i < s.size(); ++i) dnorm += s(i) * s(i);
  // Return weight
  return (snorm / (snorm + dnorm));
}

template<size_t NA, size_t K, class Q = Quantum>
double decompose
(const QSDArray<NA, Q>& a, const IVector<K>& index, QSDArray<NA-K+1, Q>& l, QSDArray<K+1, Q>& r, int D = 0)
{
  const size_t NL = NA - K + 1;
  const size_t NR = K + 1;

  IVector<NA> a_permute;
  decompose_shape(index, a_permute);

  QSDArray<NA> a_copy;
  QSDpermute(a, a_permute, a_copy);

  DiagonalQSDArray<1> s;
  double dnorm = QSDgesvd(LeftCanonical, a_copy, s, l, r, D);
  SDdidm(s, r);

  double snorm = SDdot(s, s);

  return (snorm / (snorm + dnorm));
}

template<size_t NA, size_t NL, class Q = Quantum>
double indexed_decompose
(const QSDArray<NA, Q>&      a, const IVector<NA>&      a_symbols,
       QSDArray<NL, Q>&      l, const IVector<NL>&      l_symbols,
       QSDArray<NA-NL+2, Q>& r, const IVector<NA-NL+2>& r_symbols, int D = 0)
{
  const size_t NR = NA - NL + 2;

  IVector<NA> a_permute;
  IVector<NL> l_permute;
  IVector<NR> r_permute;
  indexed_decompose_shape(a_symbols, l_symbols, r_symbols,
                          a_permute, l_permute, r_permute);

  QSDArray<NA> a_copy;
  QSDpermute(a, a_permute, a_copy);

  DiagonalQSDArray<1> s;
  QSDArray<NL> l_tmp;
  QSDArray<NL> r_tmp;
  double dnorm = QSDgesvd(LeftCanonical, a_copy, s, l_tmp, r_tmp, D);
  SDdidm(s, r_tmp);

  double snorm = SDdot(s, s);
  QSDpermute(l_tmp, l_permute, l);
  QSDpermute(r_tmp, r_permute, r);

  return (snorm / (snorm + dnorm));
}

};

#endif // _BTAS_CXX11_DRIVER_DECOMPOSE_H
