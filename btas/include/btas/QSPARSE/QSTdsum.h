#ifndef _BTAS_CXX11_QSTDSUM_H
#define _BTAS_CXX11_QSTDSUM_H 1

#include <btas/TVector.h>

#include <btas/QSPARSE/QSTArray.h>
#include <btas/SPARSE/STdsum.h>

namespace btas {

//! Direct sum of arrays X and Y, adding to Z
/*! E.g.) Z = X (+) Y : N = 2
 *  | x x |                 | x x 0 0 0 |
 *  | x x | (+)           = | x x 0 0 0 |
 *              | y y y |   | 0 0 y y y |
 */
template<typename T, size_t N, class Q = Quantum>
void QSTdsum(const QSTArray<T, N, Q>& x, const QSTArray<T, N, Q>& y, QSTArray<T, N, Q>& z)
{
  if(x.q() != y.q())
      BTAS_THROW(false, "btas::QSTdsum: total quantum numbers mismatched");
  const TVector<Qshapes<Q>, N>& x_qshape = x.qshape();
  const TVector<Qshapes<Q>, N>& y_qshape = y.qshape();
        TVector<Qshapes<Q>, N>  z_qshape;
  for(int i = 0; i < N; ++i) z_qshape[i] = x_qshape[i] + y_qshape[i];
  // Check sparse shape of z
  if(z.size() > 0) {
    if(x.q() != z.q())
      BTAS_THROW(false, "btas::QSTdsum: total quantum number of z mismatched");
    if(z.qshape() != z_qshape)
      BTAS_THROW(false, "btas::QSTdsum: quantum number indices of z mismatched");
  }
  else {
    TVector<Dshapes, N> z_dshape = x.dshape();
    for(int i = 0; i < N; ++i) z_dshape[i].insert(z_dshape[i].end(), y.dshape(i).begin(), y.dshape(i).end());
    z.resize(x.q(), z_qshape, z_dshape, false);
  }
  // Call sparse-array function
  STdsum<T, N>(x, y, z);
}

//! Partial direct sum of arrays X and Y, adding to Z
/*! E.g.) Z = X (+) Y with trace index = { 0 } : N = 2, K = 1
 *  | x x |     | y y y |   | x x y y y |
 *  | x x | (+) | y y y | = | x x y y y |
 *  | x x |     | y y y |   | x x y y y |
 *
 *  Note that sizes of trace indices must be the same between two arrays
 */
template<typename T, size_t N, size_t K, class Q = Quantum>
void QSTdsum(const QSTArray<T, N, Q>& x, const QSTArray<T, N, Q>& y, const IVector<K>& trace_index, QSTArray<T, N, Q>& z)
{
  if(x.q() != y.q())
      BTAS_THROW(false, "btas::QSTdsum: total quantum numbers mismatched");
  const TVector<Qshapes<Q>, N>& x_qshape = x.qshape();
  const TVector<Qshapes<Q>, N>& y_qshape = y.qshape();
  for(int k = 0; k < K; ++k) {
    int tk = trace_index[k];
    if(x_qshape[tk] != y_qshape[tk])
      BTAS_THROW(false, "btas::QSTdsum: found mismatched quantum number indices to be traced");
  }
  IVector<N-K> dsum_index; // Complement of trace_index
  int nsum = 0;
  for(int i = 0; i < N; ++i) {
    if(std::find(trace_index.begin(), trace_index.end(), i) == trace_index.end()) dsum_index[nsum++] = i;
  }
  TVector<Qshapes<Q>, N> z_qshape(x_qshape);
  for(int i = 0; i < N-K; ++i) z_qshape[dsum_index[i]] += y_qshape[dsum_index[i]];
  // Check sparse shape of z
  if(z.size() > 0) {
    if(x.q() != z.q())
      BTAS_THROW(false, "btas::QSTdsum: total quantum number of z mismatched");
    if(z.qshape() != z_qshape)
      BTAS_THROW(false, "btas::QSTdsum: array shape of z mismatched");
  }
  else {
    const TVector<Dshapes, N>& x_dshape = x.dshape();
    const TVector<Dshapes, N>& y_dshape = y.dshape();
          TVector<Dshapes, N>  z_dshape = x_dshape;
    for(int i = 0; i < N-K; ++i)
      z_dshape[dsum_index[i]].insert(z_dshape[dsum_index[i]].end(), y_dshape[dsum_index[i]].begin(), y_dshape[dsum_index[i]].end());
    z.resize(x.q(), z_qshape, z_dshape, false);
  }
  // Call sparse-array function
  STdsum<T, N, K>(x, y, trace_index, z);
}

}; // namespace btas

#endif // _BTAS_CXX11_QSTDSUM_H
