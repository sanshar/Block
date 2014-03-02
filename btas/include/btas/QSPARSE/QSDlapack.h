#ifndef _BTAS_CXX11_QSDLAPACK_H
#define _BTAS_CXX11_QSDLAPACK_H 1

#include <btas/btas.h>

#include <btas/DENSE/Dlapack.h>
#include <btas/DENSE/Darglist.h>

#include <btas/SPARSE/SDblas.h>

#include <btas/QSPARSE/QSDArray.h>
#include <btas/QSPARSE/QSDmerge.h>

namespace btas {

/*! \brief Arrow direction of decomposed array
 *  Because SVD splits A = U * S * V^T, either U or V has a total quantum number equals to that of A.
 *
 *  \par LeftArrow:
 *  Defined U has the quantum number of A and V^T has zero quantum number
 *
 *  \par RightArrow:
 *  Defined V^T has the quantum number of A and U has zero quantum number
 */
enum BTAS_ARROW_DIRECTION {
  LeftArrow, //!< U has the quantum number of A in SVD
  RightArrow //!< V^T has the quantum number of A in SVD
};

/*! \brief Calling Dgesvd threaded for each non-zero block
 *  Only for matrix form.
 *  Suppose matrix A is block-diagonal, i.e. quantum number indices must be merged, which simplifies computation.
 *
 *  QSDgesvd subroutine merges array A to pass it here.
 */
template<class Q = Quantum>
void thread_QSDgesvd
(const BTAS_ARROW_DIRECTION& ArrowDir,
 const QSDArray<2, Q>& a, SDArray<1>& s, QSDArray<2, Q>& u, QSDArray<2, Q>& vt, bool calc_full_u = false, bool calc_full_vt = false)
{
  int nrows = a.shape(0);
  int ncols = a.shape(1);
  // Contraction list for thread parallelism
  std::vector<DgesvdArglist<2, 2>> task_list;
  task_list.reserve(a.size());
  // Set task list
  for(typename QSDArray<2, Q>::const_iterator ia = a.begin(); ia != a.end(); ++ia) {
    if(ia->second->size() == 0) continue;
    // Calc. block indices of U, S, V^T
    int irow = ia->first / ncols;
    int icol = ia->first % ncols;
    int stag, utag, vtag;
    if(ArrowDir == LeftArrow) {
      stag = irow;
      utag = irow * nrows + irow;
      vtag = irow * ncols + icol;
    }
    else {
      stag = icol;
      utag = irow * ncols + icol;
      vtag = icol * ncols + icol;
    }
    SDArray<1>::iterator is = s.reserve(stag);
    SDArray<2>::iterator iu = u.reserve(utag);
    SDArray<2>::iterator iv = vt.reserve(vtag);
    DgesvdArglist<2, 2> gesvd_list(ia->second, is->second, iu->second, iv->second, calc_full_u, calc_full_vt);
    task_list.push_back(gesvd_list);
  }
  parallel_call(task_list);
}

/*!
 * \brief Thin Singular Value Decomposition
 *
 * If Dim = 0, all non-zero singular values (>= 1.0e-16) are kept
 * If Dim > 0, only Dim number of singular values are kept
 * If Dim < 0, discards singular values less than Tol x 10^(Dim)
 *
 * Returns total discarded norm (or density weights): sum_{i > Dim} S(i)^2
 */
template<size_t NA, size_t NU, class Q = Quantum>
double QSDgesvd
(const BTAS_ARROW_DIRECTION& ArrowDir,
 const QSDArray<NA, Q>& a, SDArray<1>& s, QSDArray<NU, Q>& u, QSDArray<NA-NU+2, Q>& vt, int Dim = 0, double Tol = 1.0)
{
  const size_t NL = NU-1;
  const size_t NR = NA-NL;

  const TVector<Qshapes<Q>, NA>& a_qshape = a.qshape();
  const TVector<Dshapes,    NA>& a_dshape = a.dshape();

  // Calc. row (left) shapes
  TVector<Qshapes<Q>, NL> a_qshape_left;
  TVector<Dshapes,    NL> a_dshape_left;

  for(int i = 0; i < NL; ++i) {
    a_qshape_left [i] = a_qshape[i];
    a_dshape_left [i] = a_dshape[i];
  }

  // Calc. col (right) shapes
  TVector<Qshapes<Q>, NR> a_qshape_right;
  TVector<Dshapes,    NR> a_dshape_right;

  for(int i = 0; i < NR; ++i) {
    a_qshape_right[i] = a_qshape[i+NL];
    a_dshape_right[i] = a_dshape[i+NL];
  }

  // Merge array A into matrix form
  QSTmergeInfo<NL, Q> a_qinfo_left (a_qshape_left,  a_dshape_left );
  QSTmergeInfo<NR, Q> a_qinfo_right(a_qshape_right, a_dshape_right);

  QSDArray<2, Q>  a_merge;
  QSDmerge(a_qinfo_left, a, a_qinfo_right, a_merge);

  // Determine arrow direction
  Qshapes<Q> q_rows(a_qinfo_left.qshape_merged());
  Qshapes<Q> q_cols(a_qinfo_right.qshape_merged());

  Q  u_q_total;
  Q vt_q_total;
  Qshapes<Q> q_sval;
  if(ArrowDir == LeftArrow) {
     u_q_total = Q::zero();
    vt_q_total = a.q();
       q_sval  = q_rows;
  }
  else
  {
     u_q_total = a.q();
    vt_q_total = Q::zero();
       q_sval  =-q_cols;
  }

  // Carry out SVD on merged matrix A
  SDArray<1>    s_value (shape(q_sval.size()));
  QSDArray<2, Q> u_merge ( u_q_total, make_array( q_rows,-q_sval));
  QSDArray<2, Q> vt_merge(vt_q_total, make_array( q_sval, q_cols));

  thread_QSDgesvd(ArrowDir, a_merge, s_value, u_merge, vt_merge);
  a_merge.clear();


  // Truncate by singular values
  // Dicarded norm: dnorm = sum_{i > D} s_value[i]^2
  double dnorm = 0.0;
  int n_sval = q_sval.size();
  int n_cols = q_cols.size();
  // Containers of selected quantum number indices and sizes
  Qshapes<Q> q_sval_nz; q_sval_nz.reserve(n_sval);
  Dshapes    d_sval_nz; d_sval_nz.reserve(n_sval);
  std::map<int, int> map_sval_nz;
  // Collect singular values
  std::vector<double> s_sorted;
  for(typename SDArray<1>::iterator its = s_value.begin(); its != s_value.end(); ++its)
    s_sorted.insert(s_sorted.end(), its->second->begin(), its->second->end());
  // Sort descending order
  std::sort(s_sorted.rbegin(), s_sorted.rend());
  // Calc. cutoff tolerance
  double cutoff = 1.0e-16;
  if(Dim > 0 && Dim < s_sorted.size())
    cutoff = s_sorted[Dim-1];
  if(Dim < 0)
    cutoff = fabs(Tol) * pow(10.0, Dim);
  // Select singular values
  int nnz = 0;
  for(typename SDArray<1>::iterator its = s_value.begin(); its != s_value.end(); ++its) {
    typename DArray<1>::iterator itd = its->second->begin();
    int D_kept = 0;
    for(; itd != its->second->end(); ++itd) {
      if(*itd < cutoff) break;
      ++D_kept;
    }
    for(; itd != its->second->end(); ++itd) {
      dnorm += (*itd) * (*itd);
    }
    if(D_kept > 0) {
      q_sval_nz.push_back(q_sval[its->first]);
      d_sval_nz.push_back(D_kept);
      map_sval_nz.insert(std::make_pair(its->first, nnz++));
    }
  }
  // Copy selected singular values
  SDArray<1> s_value_nz(shape(nnz));
  for(typename SDArray<1>::iterator it = s_value.begin(); it != s_value.end(); ++it) {
    typename std::map<int, int>::iterator imap = map_sval_nz.find(it->first);
    if(imap != map_sval_nz.end()) {
      typename SDArray<1>::iterator jt = s_value_nz.reserve(imap->second);
      int Ds = d_sval_nz[imap->second];
      jt->second->resize(Ds);
      *jt->second = it->second->subarray(shape(0), shape(Ds-1));
    }
  }
  s_value_nz.check_dshape();
  s_value.clear();
  // Copy selected left-singular vectors
  QSDArray<2, Q> u_merge_nz(u_merge.q(), make_array( q_rows,-q_sval_nz));
  for(typename QSDArray<2, Q>::iterator it = u_merge.begin(); it != u_merge.end(); ++it) {
    int irow = it->first / n_sval;
    int icol = it->first % n_sval;
    typename std::map<int, int>::iterator imap = map_sval_nz.find(icol);
    if(imap != map_sval_nz.end()) {
      typename QSDArray<2, Q>::iterator jt = u_merge_nz.reserve(irow * nnz + imap->second);
      assert(jt != u_merge_nz.end()); // if aborted here, there's a bug in btas::QSDgesvd
      int Ds = d_sval_nz[imap->second];
      int Dr = it->second->shape(0);
      jt->second->resize(Dr, Ds);
      *jt->second = it->second->subarray(shape(0, 0), shape(Dr-1, Ds-1));
    }
  }
  u_merge_nz.check_dshape();
  u_merge.clear();
  // Copy selected right-singular vectors
  QSDArray<2, Q> vt_merge_nz(vt_merge.q(), make_array( q_sval_nz, q_cols));
  for(typename QSDArray<2, Q>::iterator it = vt_merge.begin(); it != vt_merge.end(); ++it) {
    int irow = it->first / n_cols;
    int icol = it->first % n_cols;
    typename std::map<int, int>::iterator imap = map_sval_nz.find(irow);
    if(imap != map_sval_nz.end()) {
      typename QSDArray<2, Q>::iterator jt = vt_merge_nz.reserve(imap->second * n_cols + icol);
      assert(jt != vt_merge_nz.end()); // if aborted here, there's a bug in btas::QSDgesvd
      int Ds = d_sval_nz[imap->second];
      int Dc = it->second->shape(1);
      jt->second->resize(Ds, Dc);
      *jt->second = it->second->subarray(shape(0, 0), shape(Ds-1, Dc-1));
    }
  }
  vt_merge_nz.check_dshape();
  vt_merge.clear();

  // Reshape matrix to array form
  SDcopy  (s_value_nz, s);
  QSDexpand(a_qinfo_left, u_merge_nz, u);
  QSDexpand(vt_merge_nz, a_qinfo_right, vt);

  return dnorm;

}

//! Full SVD
/*!
 *  \param s_rm removed singular values
 *  \param u_rm null space of left singular vectors
 *  \param vt_rm null space of right singular vectors
 */
template<size_t NA, size_t NU, class Q = Quantum>
   double QSDgesvd
(const BTAS_ARROW_DIRECTION& ArrowDir,
 const QSDArray<NA, Q>& a, SDArray<1>& s_nz, SDArray<1>& s_rm,
 bool calc_full_u, QSDArray<NU, Q>& u_nz, QSDArray<NU, Q>& u_rm,
 bool calc_full_vt, QSDArray<NA-NU+2, Q>& vt_nz, QSDArray<NA-NU+2, Q>& vt_rm, int Dim = 0, double Tol = 1.0)
{
   const size_t NL = NU-1;
   const size_t NR = NA-NL;
   const TVector<Qshapes<Q>, NA>& a_qshape = a.qshape();
   const TVector<Dshapes,    NA>& a_dshape = a.dshape();
   // Calc. row (left) shapes
   TVector<Qshapes<Q>, NL> a_qshape_left;
   TVector<Dshapes,    NL> a_dshape_left;
   for(int i = 0; i < NL; ++i) {
      a_qshape_left [i] = a_qshape[i];
      a_dshape_left [i] = a_dshape[i];
   }
   // Calc. col (right) shapes
   TVector<Qshapes<Q>, NR> a_qshape_right;
   TVector<Dshapes,    NR> a_dshape_right;
   for(int i = 0; i < NR; ++i) {
      a_qshape_right[i] = a_qshape[i+NL];
      a_dshape_right[i] = a_dshape[i+NL];
   }
   // Merge array A into matrix form
   QSTmergeInfo<NL, Q> a_qinfo_left (a_qshape_left,  a_dshape_left );
   QSTmergeInfo<NR, Q> a_qinfo_right(a_qshape_right, a_dshape_right);
   QSDArray<2, Q>  a_merge;
   QSDmerge(a_qinfo_left, a, a_qinfo_right, a_merge);
   // Determine arrow direction
   Qshapes<Q> q_rows(a_qinfo_left.qshape_merged());
   Qshapes<Q> q_cols(a_qinfo_right.qshape_merged());
   Q  u_q_total;
   Q vt_q_total;
   Qshapes<Q> q_sval;
   if(ArrowDir == LeftArrow) {
      u_q_total = Q::zero();
      vt_q_total = a.q();
      q_sval  = q_rows;
   }
   else
   {
      u_q_total = a.q();
      vt_q_total = Q::zero();
      q_sval  =-q_cols;
   }
   // Carry out SVD on merged matrix A
   SDArray<1>    s_value (shape(q_sval.size()));
   QSDArray<2, Q> u_merge ( u_q_total, make_array( q_rows,-q_sval));
   QSDArray<2, Q> vt_merge(vt_q_total, make_array( q_sval, q_cols));
   thread_QSDgesvd(ArrowDir, a_merge, s_value, u_merge, vt_merge, calc_full_u, calc_full_vt);
   a_merge.clear();

   // Truncate by singular values
   // Dicarded norm: dnorm = sum_{i > D} s_value[i]^2
   double dnorm = 0.0;
   int n_sval = q_sval.size();
   int n_cols = q_cols.size();
   // Containers of selected quantum number indices and sizes
   Qshapes<Q> q_sval_nz; q_sval_nz.reserve(n_sval);
   Dshapes    d_sval_nz; d_sval_nz.reserve(n_sval);
   std::map<int, int> map_sval_nz;
   // Containers of removed quantum number indices and sizes
   Qshapes<Q> q_sval_rm; q_sval_rm.reserve(n_sval);
   Dshapes    d_sval_rm; d_sval_rm.reserve(n_sval);
   std::map<int, int> map_sval_rm;
   // Collect singular values
   std::vector<double> s_sorted;
   for(typename SDArray<1>::iterator its = s_value.begin(); its != s_value.end(); ++its)
      s_sorted.insert(s_sorted.end(), its->second->begin(), its->second->end());
   // Sort descending order
   std::sort(s_sorted.rbegin(), s_sorted.rend());
   // Calc. cutoff tolerance
   double cutoff = 1.0e-16;
   if(Dim > 0 && Dim < s_sorted.size())
      cutoff = s_sorted[Dim-1];
   if(Dim < 0)
      cutoff = fabs(Tol) * pow(10.0, Dim);
   // Select singular values
   int nnz = 0;
   int nrm = 0;
   for(typename SDArray<1>::iterator its = s_value.begin(); its != s_value.end(); ++its) {
      typename DArray<1>::iterator itd = its->second->begin();
      int D_kept = 0;
      for(; itd != its->second->end(); ++itd) {
         if(*itd < cutoff) break;
         ++D_kept;
      }
      int D_remv = its->second->size()-D_kept;

      if(D_kept > 0) {
         q_sval_nz.push_back(q_sval[its->first]);
         d_sval_nz.push_back(D_kept);
         map_sval_nz.insert(std::make_pair(its->first, nnz++));
      }
      if(D_remv > 0) {
         q_sval_rm.push_back(q_sval[its->first]);
         d_sval_rm.push_back(D_remv);
         map_sval_rm.insert(std::make_pair(its->first, nrm++));
      }
   }
   // Copying singular values
   SDArray<1> s_value_nz(shape(nnz));
   SDArray<1> s_value_rm(shape(nrm));
   for(typename SDArray<1>::iterator it = s_value.begin(); it != s_value.end(); ++it) {
      int Ds = 0;
      typename std::map<int, int>::iterator imap = map_sval_nz.find(it->first);
      if(imap != map_sval_nz.end()) {
         typename SDArray<1>::iterator jt = s_value_nz.reserve(imap->second);
         Ds = d_sval_nz[imap->second];
         jt->second->resize(Ds);
         *jt->second = it->second->subarray(shape(0), shape(Ds-1));
      }
      typename std::map<int, int>::iterator jmap = map_sval_rm.find(it->first);
      if(jmap != map_sval_rm.end()) {
         typename SDArray<1>::iterator jt = s_value_rm.reserve(jmap->second);
         int Dx = d_sval_rm[jmap->second];
         jt->second->resize(Dx);
         *jt->second = it->second->subarray(shape(Ds), shape(Ds+Dx-1));
      }
   }
   s_value_nz.check_dshape();
   s_value_rm.check_dshape();
   s_value.clear();
   // Copying left-singular vectors
   QSDArray<2, Q> u_merge_nz(u_merge.q(), make_array( q_rows,-q_sval_nz));
   QSDArray<2, Q> u_merge_rm(u_merge.q(), make_array( q_rows,-q_sval_rm));
   for(typename QSDArray<2, Q>::iterator it = u_merge.begin(); it != u_merge.end(); ++it) {
      int irow = it->first / n_sval;
      int icol = it->first % n_sval;
      int Ds = 0;
      int Dr = it->second->shape(0);
      typename std::map<int, int>::iterator imap = map_sval_nz.find(icol);
      if(imap != map_sval_nz.end()) {
         typename QSDArray<2, Q>::iterator jt = u_merge_nz.reserve(irow * nnz + imap->second);
         assert(jt != u_merge_nz.end()); // if aborted here, there's a bug in btas::QSDgesvd
         Ds = d_sval_nz[imap->second];
         jt->second->resize(Dr, Ds);
         *jt->second = it->second->subarray(shape(0, 0), shape(Dr-1, Ds-1));
      }
      typename std::map<int, int>::iterator jmap = map_sval_rm.find(icol);
      if(jmap != map_sval_rm.end()) {
         typename QSDArray<2, Q>::iterator jt = u_merge_rm.reserve(irow * nrm + jmap->second);
         assert(jt != u_merge_rm.end()); // if aborted here, there's a bug in btas::QSDgesvd
         int Dx = d_sval_rm[jmap->second];
         jt->second->resize(Dr, Dx);
         *jt->second = it->second->subarray(shape(0, Ds), shape(Dr-1, Ds+Dx-1));
      }
   }
   u_merge_nz.check_dshape();
   u_merge_rm.check_dshape();
   u_merge.clear();
   // Copy selected right-singular vectors
   QSDArray<2, Q> vt_merge_nz(vt_merge.q(), make_array( q_sval_nz, q_cols));
   QSDArray<2, Q> vt_merge_rm(vt_merge.q(), make_array( q_sval_rm, q_cols));
   for(typename QSDArray<2, Q>::iterator it = vt_merge.begin(); it != vt_merge.end(); ++it) {
      int irow = it->first / n_cols;
      int icol = it->first % n_cols;
      int Ds = 0;
      int Dc = it->second->shape(1);
      typename std::map<int, int>::iterator imap = map_sval_nz.find(irow);
      if(imap != map_sval_nz.end()) {
         typename QSDArray<2, Q>::iterator jt = vt_merge_nz.reserve(imap->second * n_cols + icol);
         assert(jt != vt_merge_nz.end()); // if aborted here, there's a bug in btas::QSDgesvd
         Ds = d_sval_nz[imap->second];
         jt->second->resize(Ds, Dc);
         *jt->second = it->second->subarray(shape(0, 0), shape(Ds-1, Dc-1));
      }
      typename std::map<int, int>::iterator jmap = map_sval_rm.find(irow);
      if(jmap != map_sval_rm.end()) {
         typename QSDArray<2, Q>::iterator jt = vt_merge_rm.reserve(jmap->second * n_cols + icol);
         assert(jt != vt_merge_rm.end()); // if aborted here, there's a bug in btas::QSDgesvd
         int Dx = d_sval_rm[jmap->second];
         jt->second->resize(Dx, Dc);
         *jt->second = it->second->subarray(shape(Ds, 0), shape(Ds+Dx-1, Dc-1));
      }
   }
   vt_merge_nz.check_dshape();
   vt_merge_rm.check_dshape();
   vt_merge.clear();
   // Reshape matrix to array form
   if(nnz > 0) {
      SDcopy  (s_value_nz, s_nz);
      QSDexpand(a_qinfo_left, u_merge_nz, u_nz);
      QSDexpand(vt_merge_nz, a_qinfo_right, vt_nz);
   }
   if(nrm > 0) {
      SDcopy  (s_value_rm, s_rm);
      QSDexpand(a_qinfo_left, u_merge_rm, u_rm);
      QSDexpand(vt_merge_rm, a_qinfo_right, vt_rm);
   }

   return SDdot(s_rm, s_rm);
}

}; // namespace btas

#endif // _BTAS_CXX11_QSDLAPACK_H
