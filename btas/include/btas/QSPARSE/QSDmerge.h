
//
/*! \file  QSDmerge.h
 *  \brief Aliases of QSTmerge and QSTexpand for double precision real array
 */

#ifndef _BTAS_CXX11_QSDMERGE_H
#define _BTAS_CXX11_QSDMERGE_H 1

#include <btas/QSPARSE/QSDArray.h>
#include <btas/QSPARSE/QSTmerge.h>

namespace btas {

template<size_t MR, size_t N, class Q = Quantum>
inline void QSDmerge
(const QSTmergeInfo<MR, Q>& rows_info, const QSDArray<N, Q>& a, QSDArray<1+N-MR, Q>& b)
{
  QSTmerge<double, MR, N, Q>(rows_info, a, b);
}

template<size_t N, size_t MC, class Q = Quantum>
inline void QSDmerge
(const QSDArray<N, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSDArray<N-MC+1, Q>& b)
{
  QSTmerge<double, N, MC, Q>(a, cols_info, b);
}

template<size_t MR, size_t MC, class Q = Quantum>
inline void QSDmerge
(const QSTmergeInfo<MR, Q>& rows_info, const QSDArray<MR+MC, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSDArray<2, Q>& b)
{
  QSTmerge<double, MR, MC, Q>(rows_info, a, cols_info, b);
}

template<size_t MR, size_t N, class Q = Quantum>
inline void QSDexpand
(const QSTmergeInfo<MR, Q>& rows_info, const QSDArray<1+N-MR, Q>& a, QSDArray<N, Q>& b)
{
  QSTexpand<double, MR, N, Q>(rows_info, a, b);
}

template<size_t N, size_t MC, class Q = Quantum>
inline void QSDexpand
(const QSDArray<N-MC+1, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSDArray<N, Q>& b)
{
  QSTexpand<double, N, MC, Q>(a, cols_info, b);
}

template<size_t MR, size_t MC, class Q = Quantum>
inline void QSDexpand
(const QSTmergeInfo<MR, Q>& rows_info, const QSDArray<2, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSDArray<MR+MC, Q>& b)
{
  QSTexpand<double, MR, MC, Q>(rows_info, a, cols_info, b);
}

};

#endif // _BTAS_CXX11_QSDMERGE_H
