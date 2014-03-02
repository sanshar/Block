#ifndef _BTAS_CXX11_SDBLAS_H
#define _BTAS_CXX11_SDBLAS_H 1

#include <algorithm>
#include <functional>
#include <cmath>


#include <btas/DENSE/Dblas.h>
#include <btas/DENSE/Darglist.h>

#include <btas/SPARSE/SDArray.h>
#include <btas/SPARSE/btas_contract_dshape.h>

#ifndef _SERIAL_REPLICATION_LIMIT
#define _SERIAL_REPLICATION_LIMIT 1
#endif

#ifndef _SERIAL_CONTRACTION_LIMIT
#define _SERIAL_CONTRACTION_LIMIT 1
#endif

namespace btas {

//####################################################################################################
// Sparse BLAS wrappers for double-precision real array : serial version
//####################################################################################################

//! Call Dcopy per every non-zero block
/*! if _up_cast is on,  only blocks in 'x' that allowed in 'y' are copied
 *  in other words, e.g. doing up-cast from SDArray(x) to QSDArray(y)
 *
 *  if _up_cast is off, it's aborted when block in 'x' that is not allowed in 'y' is found
 */
template<size_t N>
void serial_SDcopy
(const SDArray<N>& x, SDArray<N>& y, bool _up_cast = false)
{
  for(typename SDArray<N>::const_iterator ix = x.begin(); ix != x.end(); ++ix) {
    typename SDArray<N>::iterator iy = y.reserve(ix->first);
    if(iy != y.end())
      Dcopy(*(ix->second), *(iy->second));
    else
      BTAS_THROW(_up_cast, "btas::serial_SDcopy; requested block must be zero, unable to be reserved");
  }
}

//! Call Ddot per every non-zero block
template<size_t N>
double serial_SDdot
(const SDArray<N>& x, const SDArray<N>& y)
{
  double dotxy = 0.0;
  for(typename SDArray<N>::const_iterator ix = x.begin(); ix != x.end(); ++ix) {
    typename SDArray<N>::const_iterator iy = y.find(ix->first);
    if(iy == y.end()) continue;
    dotxy += Ddot(*(ix->second), *(iy->second));
  }
  return dotxy;
}

//! Call Dscal per every non-zero block
template<size_t N>
void serial_SDscal
(const double& alpha, SDArray<N>& x)
{
  for(typename SDArray<N>::iterator ix = x.begin(); ix != x.end(); ++ix)
    Dscal(alpha, *(ix->second));
}

//! Call Daxpy per every non-zero block
template<size_t N>
void serial_SDaxpy
(const double& alpha, const SDArray<N>& x, SDArray<N>& y) {
  for(typename SDArray<N>::const_iterator ix = x.begin(); ix != x.end(); ++ix)
{
    typename SDArray<N>::iterator iy = y.reserve(ix->first);
    if(iy == y.end())
      BTAS_THROW(false, "btas::serial_SDaxpy; requested block must be zero, unable to be reserved");
    Daxpy(alpha, *(ix->second), *(iy->second));
  }
}

//####################################################################################################
// Sparse BLAS wrappers for double-precision real array : threaded by OpenMP
//####################################################################################################

//! Call Dcopy per every non-zero block
/*! if _up_cast is on,  only blocks in 'x' that allowed in 'y' are copied
 *  in other words, e.g. doing up-cast from SDArray(x) to QSDArray(y)
 *
 *  if _up_cast is off, it's aborted when block in 'x' that is not allowed in 'y' is found
 */
template<size_t N>
void thread_SDcopy
(const SDArray<N>& x, SDArray<N>& y, bool _up_cast)
{
  std::vector<DcopyArglist<N>> task_list;
  task_list.reserve(x.nnz());
  for(typename SDArray<N>::const_iterator ix = x.begin(); ix != x.end(); ++ix) {
    typename SDArray<N>::iterator iy = y.reserve(ix->first);
    if(iy != y.end())
      task_list.push_back(DcopyArglist<N>(ix->second, iy->second));
    else
      BTAS_THROW(_up_cast, "btas::thread_SDcopy; requested block must be zero, unable to be reserved");
  }
  parallel_call(task_list);
}

//! Call Dscal per every non-zero block
template<size_t N>
void thread_SDscal
(const double& alpha, SDArray<N>& x)
{
  std::vector<DscalArglist<N>> task_list;
  task_list.reserve(x.nnz());
  for(typename SDArray<N>::iterator ix = x.begin(); ix != x.end(); ++ix)
    task_list.push_back(DscalArglist<N>(alpha, ix->second));
  parallel_call(task_list);
}

//! Call Daxpy per every non-zero block
template<size_t N>
void thread_SDaxpy
(const double& alpha, const SDArray<N>& x, SDArray<N>& y)
{
  std::vector<DaxpyArglist<N>> task_list;
  task_list.reserve(x.nnz());
  for(typename SDArray<N>::const_iterator ix = x.begin(); ix != x.end(); ++ix) {
    typename SDArray<N>::iterator iy = y.reserve(ix->first);
    if(iy == y.end())
      BTAS_THROW(false, "btas::thread_SDaxpy; requested block must be zero, unable to be reserved");
    task_list.push_back(DaxpyArglist<N>(alpha, ix->second, iy->second));
  }
  parallel_call(task_list);
}

//! Call Dgemv per every non-zero block
/*! before calling thread_SDgemv 
 *  	- 'c' is scaled by beta
 *  	- sparse-block of 'a' is transposed if TransA = Trans,
 *        note that dense-arrays have not yet been transposed
 */
template<size_t NA, size_t NB, size_t NC>
void thread_SDgemv
(const BTAS_TRANSPOSE& TransA,
 const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, SDArray<NC>& c)
{
  int nrows  = std::accumulate(a.shape().begin(), a.shape().begin()+NC, 1, std::multiplies<int>());
  int stride = std::accumulate(b.shape().begin(), b.shape().end(),      1, std::multiplies<int>());
  // contraction list for thread parallelism
  std::vector<DgemvArglist<NA, NB, NC>> task_list;
  task_list.reserve(a.nnz());
  // block contraction
  for(int i = 0; i < nrows; ++i) {
    // calc. range in 'a' contributes to c(i)
    int a_lbound = i * stride;
    int a_ubound = a_lbound + stride - 1;
    typename SDArray<NA>::const_iterator ialo = a.lower_bound(a_lbound);
    typename SDArray<NA>::const_iterator iaup = a.upper_bound(a_ubound);
    if(ialo == iaup)  continue;
    if(!c.allowed(i)) continue;
    // calc. find block in 'a' & 'b' contributes to c(i) to make argment list
    DgemvArglist<NA, NB, NC> gemv_args;
    for(typename SDArray<NA>::const_iterator ia = ialo; ia != iaup; ++ia) {
      typename SDArray<NB>::const_iterator jb = b.find(ia->first % stride);
      if(jb != b.end())
        gemv_args.add(ia->second, jb->second);
    }
    // there's no blocks in 'a', 'b' contribute to c(i)
    if(gemv_args.size() == 0) continue;
    // allocate block element for c(i)
    typename SDArray<NC>::iterator ic = c.reserve(i);
    if(ic == c.end())
      BTAS_THROW(false, "btas::thread_SDgemv required block could not be allocated");
    // add argments to task_list
    gemv_args.reset(ic->second, TransA, alpha, 1.0);
    task_list.push_back(gemv_args);
  }
  parallel_call(task_list);
}

//! Call Dger per every non-zero block
template<size_t NA, size_t NB, size_t NC>
void thread_SDger
(const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, SDArray<NC>& c)
{
  int stride = std::accumulate(b.shape().begin(), b.shape().end(), 1, std::multiplies<int>());
  // contraction list for thread parallelism
  std::vector<DgerArglist<NA, NB, NC>> task_list;
  task_list.reserve(a.nnz()*b.nnz()); /*! approximately */
  // block contraction
  for(typename SDArray<NA>::const_iterator ia = a.begin(); ia != a.end(); ++ia) {
    int c_irow = ia->first * stride;
    for(typename SDArray<NB>::const_iterator jb = b.begin(); jb != b.end(); ++jb) {
      int c_tag = c_irow + jb->first;
      if(!c.allowed(c_tag)) continue;
      // allocate block element for c(i, j)
      typename SDArray<NC>::iterator ic = c.reserve(c_tag);
      if(ic == c.end())
        BTAS_THROW(false, "btas::thread_SDger required block could not be allocated");
      // make argment list
      DgerArglist<NA, NB, NC> ger_args(ic->second, alpha);
      ger_args.add(ia->second, jb->second);
      task_list.push_back(ger_args);
    }
  }
  parallel_call(task_list);
}

//! Call Dgemm per every non-zero block
/*! before calling thread_SDgemv 
 *  	- 'c' is scaled by beta
 *  	- sparse-block of 'a' is transposed if TransA = Trans
 *      - sparse-block of 'b' is transposed if TransB = NoTrans
 *  since c(i, j) = sum_{k} a(i, k) * b(k, j), it's advantageous to store 'b' as b(j, k) order
 *  so that summation over k has less memory-miss when searching non-zero element
 *
 *  note that dense-arrays have not yet been transposed
 */
template<size_t NA, size_t NB, size_t NC>
void thread_SDgemm
(const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
 const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, SDArray<NC>& c)
{
  const size_t K = (NA + NB - NC)/2;
  int nrows  = std::accumulate(a.shape().begin(), a.shape().begin()+NA-K, 1, std::multiplies<int>());
  int stride = std::accumulate(a.shape().begin()+NA-K, a.shape().end(),   1, std::multiplies<int>());
  int ncols  = std::accumulate(b.shape().begin(), b.shape().begin()+NB-K, 1, std::multiplies<int>());
  // contraction list for thread parallelism
  std::vector<DgemmArglist<NA, NB, NC>> task_list;
  task_list.reserve(std::max(a.nnz(), b.nnz())); /*! approximately */
  // block contraction
  for(int i = 0; i < nrows; ++i) {
    // calc. range in 'a' contributes to c(i, j = 0...ncols)
    int a_lbound = i * stride;
    int a_ubound = a_lbound + stride - 1;
    typename SDArray<NA>::const_iterator ialo = a.lower_bound(a_lbound);
    typename SDArray<NA>::const_iterator iaup = a.upper_bound(a_ubound);
    if(ialo == iaup) continue;
    // calc. tag for c(i, j)
    int c_irow = i * ncols;
    for(int j = 0; j < ncols; ++j) {
      int c_tag = c_irow + j;
      if(!c.allowed(c_tag)) continue;
      // calc. range in 'a' contributes to c(i, j)
      int b_lbound = j * stride;
      int b_ubound = b_lbound + stride - 1;
      typename SDArray<NB>::const_iterator jblo = b.lower_bound(b_lbound);
      typename SDArray<NB>::const_iterator jbup = b.upper_bound(b_ubound);
      if(jblo == jbup) continue;
      // make argment list
      DgemmArglist<NA, NB, NC> gemm_args;
      for(typename SDArray<NA>::const_iterator ia = ialo; ia != iaup; ++ia) {
        for(typename SDArray<NB>::const_iterator jb = jblo; jb != jbup; ++jb) {
          if((ia->first % stride) == (jb->first % stride))
            gemm_args.add(ia->second, jb->second);
        }
      }
      if(gemm_args.size() == 0) continue;
      // allocate block element for c(i, j)
      typename SDArray<NC>::iterator ic = c.reserve(c_tag);
      if(ic == c.end())
        BTAS_THROW(false, "btas::thread_SDgemm required block could not be allocated");
      // add argments to task_list
      gemm_args.reset(ic->second, TransA, TransB, alpha, 1.0);
      task_list.push_back(gemm_args);
    }
  }
  parallel_call(task_list);
}

//####################################################################################################
// Sparse BLAS wrappers for double-precision real array : index-based scaling
//####################################################################################################

//! Call Dgemv per every non-zero block
/*! before calling thread_SDgemv 
 *  	- 'c' is scaled by beta
 *  	- sparse-block of 'a' is transposed if TransA = Trans,
 *        note that dense-arrays have not yet been transposed
 */
template<size_t NA, size_t NB, size_t NC>
void thread_SDgemv
(function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)> f_scale,
 const BTAS_TRANSPOSE& TransA,
 const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, SDArray<NC>& c)
{
  int nrows  = std::accumulate(a.shape().begin(), a.shape().begin()+NC, 1, std::multiplies<int>());
  int stride = std::accumulate(b.shape().begin(), b.shape().end(),      1, std::multiplies<int>());
  // contraction list for thread parallelism
  std::vector<DgemvArglist<NA, NB, NC>> task_list;
  task_list.reserve(a.nnz());
  // block contraction
  for(int i = 0; i < nrows; ++i) {
    // calc. range in 'a' contributes to c(i)
    int a_lbound = i * stride;
    int a_ubound = a_lbound + stride - 1;
    typename SDArray<NA>::const_iterator ialo = a.lower_bound(a_lbound);
    typename SDArray<NA>::const_iterator iaup = a.upper_bound(a_ubound);
    if(ialo == iaup)  continue;
    if(!c.allowed(i)) continue;
    // calc. find block in 'a' & 'b' contributes to c(i) to make argment list
    IVector<NC> c_index = c.index(i);
    DgemvArglist<NA, NB, NC> gemv_args;
    for(typename SDArray<NA>::const_iterator ia = ialo; ia != iaup; ++ia) {
      typename SDArray<NB>::const_iterator jb = b.find(ia->first % stride);
      if(jb != b.end())
        gemv_args.add(ia->second, jb->second, f_scale(a.index(ia->first), b.index(jb->first), c_index));
    }
    // there's no blocks in 'a', 'b' contribute to c(i)
    if(gemv_args.size() == 0) continue;
    // allocate block element for c(i)
    typename SDArray<NC>::iterator ic = c.reserve(i);
    if(ic == c.end())
      BTAS_THROW(false, "btas::thread_SDgemv required block could not be allocated");
    // add argments to task_list
    gemv_args.reset(ic->second, TransA, alpha, 1.0);
    task_list.push_back(gemv_args);
  }
  parallel_call(task_list);
}

//! Call Dger per every non-zero block
template<size_t NA, size_t NB, size_t NC>
void thread_SDger
(function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)> f_scale,
 const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, SDArray<NC>& c)
{
  int stride = std::accumulate(b.shape().begin(), b.shape().end(), 1, std::multiplies<int>());
  // contraction list for thread parallelism
  std::vector<DgerArglist<NA, NB, NC>> task_list;
  task_list.reserve(a.nnz()*b.nnz()); /*! approximately */
  // block contraction
  for(typename SDArray<NA>::const_iterator ia = a.begin(); ia != a.end(); ++ia) {
    int c_irow = ia->first * stride;
    IVector<NA> a_index = a.index(ia->first);
    for(typename SDArray<NB>::const_iterator jb = b.begin(); jb != b.end(); ++jb) {
      int c_tag = c_irow + jb->first;
      if(!c.allowed(c_tag)) continue;
      // allocate block element for c(i, j)
      typename SDArray<NC>::iterator ic = c.reserve(c_tag);
      if(ic == c.end())
        BTAS_THROW(false, "btas::thread_SDger required block could not be allocated");
      // make argment list
      DgerArglist<NA, NB, NC> ger_args(ic->second, alpha);
      ger_args.add(ia->second, jb->second, f_scale(a_index, b.index(jb->first), c.index(c_tag)));
      task_list.push_back(ger_args);
    }
  }
  parallel_call(task_list);
}

//! Call Dgemm per every non-zero block
/*! before calling thread_SDgemv 
 *  	- 'c' is scaled by beta
 *  	- sparse-block of 'a' is transposed if TransA = Trans
 *      - sparse-block of 'b' is transposed if TransB = NoTrans
 *  since c(i, j) = sum_{k} a(i, k) * b(k, j), it's advantageous to store 'b' as b(j, k) order
 *  so that summation over k has less memory-miss when searching non-zero element
 *
 *  note that dense-arrays have not yet been transposed
 */
template<size_t NA, size_t NB, size_t NC>
void thread_SDgemm
(function<double(const IVector<NA>&, const IVector<NB>&, const IVector<NC>&)> f_scale,
 const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
 const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, SDArray<NC>& c)
{
  const size_t K = (NA + NB - NC)/2;
  int nrows  = std::accumulate(a.shape().begin(), a.shape().begin()+NA-K, 1, std::multiplies<int>());
  int stride = std::accumulate(a.shape().begin()+NA-K, a.shape().end(),   1, std::multiplies<int>());
  int ncols  = std::accumulate(b.shape().begin(), b.shape().begin()+NB-K, 1, std::multiplies<int>());
  // contraction list for thread parallelism
  std::vector<DgemmArglist<NA, NB, NC>> task_list;
  task_list.reserve(std::max(a.nnz(), b.nnz())); /*! approximately */
  // block contraction
  for(int i = 0; i < nrows; ++i) {
    // calc. range in 'a' contributes to c(i, j = 0...ncols)
    int a_lbound = i * stride;
    int a_ubound = a_lbound + stride - 1;
    typename SDArray<NA>::const_iterator ialo = a.lower_bound(a_lbound);
    typename SDArray<NA>::const_iterator iaup = a.upper_bound(a_ubound);
    if(ialo == iaup) continue;
    // calc. tag for c(i, j)
    int c_irow = i * ncols;
    for(int j = 0; j < ncols; ++j) {
      int c_tag = c_irow + j;
      if(!c.allowed(c_tag)) continue;
      IVector<NC> c_index = c.index(c_tag);
      // calc. range in 'a' contributes to c(i, j)
      int b_lbound = j * stride;
      int b_ubound = b_lbound + stride - 1;
      typename SDArray<NB>::const_iterator jblo = b.lower_bound(b_lbound);
      typename SDArray<NB>::const_iterator jbup = b.upper_bound(b_ubound);
      if(jblo == jbup) continue;
      // make argment list
      DgemmArglist<NA, NB, NC> gemm_args;
      for(typename SDArray<NA>::const_iterator ia = ialo; ia != iaup; ++ia) {
        IVector<NA> a_index = a.index(ia->first);
        for(typename SDArray<NB>::const_iterator jb = jblo; jb != jbup; ++jb) {
          if((ia->first % stride) == (jb->first % stride))
            gemm_args.add(ia->second, jb->second, f_scale(a_index, b.index(jb->first), c_index));
        }
      }
      if(gemm_args.size() == 0) continue;
      // allocate block element for c(i, j)
      typename SDArray<NC>::iterator ic = c.reserve(c_tag);
      if(ic == c.end())
        BTAS_THROW(false, "btas::thread_SDgemm required block could not be allocated");
      // add argments to task_list
      gemm_args.reset(ic->second, TransA, TransB, alpha, 1.0);
      task_list.push_back(gemm_args);
    }
  }
  parallel_call(task_list);
}

//####################################################################################################
// Convenient Sparse BLAS wrapper for double precision real array
//####################################################################################################

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BLAS LEVEL 1
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<size_t N>
void SDcopy
(const SDArray<N>& x, SDArray<N>& y, bool _up_cast = false)
{
  if(_up_cast) {
    // if up-casting, 'y' must be block-shaped to be the same as 'x'
    if(x.shape() != y.shape())
      BTAS_THROW(false, "btas::SDcopy; array shape mismatched despite up-casting was specified");
  }

  y.resize(x.dshape(), false);

#ifdef _SERIAL
  serial_SDcopy(x, y, _up_cast);
#else
  thread_SDcopy(x, y, _up_cast);
#endif
}

template<size_t N>
void SDscal
(const double& alpha, SDArray<N>& x)
{
#ifdef _SERIAL
  serial_SDscal(alpha, x);
#else
  thread_SDscal(alpha, x);
#endif
}

template<size_t N>
double SDdot
(const SDArray<N>& x, const SDArray<N>& y)
{
  if(x.shape() != y.shape())
    BTAS_THROW(false, "btas::SDdot: shapes of x and y mismatched");
  return serial_SDdot(x, y);
}

template<size_t N>
void SDaxpy
(const double& alpha, const SDArray<N>& x, SDArray<N>& y)
{
  if(y.size() > 0) {
    if(x.dshape() != y.dshape())
      BTAS_THROW(false, "btas::SDaxpy: shape of y mismatched");
  }
  else {
    y.resize(x.dshape(), false);
  }
#ifdef _SERIAL
  serial_SDaxpy(alpha, x, y);
#else
  thread_SDaxpy(alpha, x, y);
#endif
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BLAS LEVEL 2
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<size_t NA, size_t NB, size_t NC>
void SDgemv
(const BTAS_TRANSPOSE& TransA,
 const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, const double& beta, SDArray<NC>& c)
{
  // calc. contraction shape
  TVector<Dshapes, NC> c_dn_shape;
  gemv_contract_dshape(TransA, a.dshape(), b.dshape(), c_dn_shape);
  if(c.size() > 0) {
    if(c_dn_shape != c.dshape())
      BTAS_THROW(false, "btas::SDgemv: block shape of c mismatched");
    SDscal(beta, c);
  }
  else {
    c.resize(c_dn_shape, false);
  }
  // call threaded sparse-dgemv
  if(TransA == NoTrans)
    thread_SDgemv(TransA, alpha, a, b, c);
  else
    thread_SDgemv(TransA, alpha, a.transposed_view(NB), b, c);
}

template<size_t NA, size_t NB, size_t NC>
void SDger
(const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, SDArray<NC>& c)
{
  // calc. contraction shape
  TVector<Dshapes, NC> c_dn_shape;
  ger_contract_dshape(a.dshape(), b.dshape(), c_dn_shape);
  if(c.size() > 0) {
    if(c_dn_shape != c.dshape())
      BTAS_THROW(false, "btas::SDger: block shape of c mismatched");
  }
  else {
    c.resize(c_dn_shape, false);
  }
  // call threaded sparse-dger
  thread_SDger(alpha, a, b, c);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BLAS level 3
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<size_t NA, size_t NB, size_t NC>
void SDgemm
(const BTAS_TRANSPOSE& TransA, const BTAS_TRANSPOSE& TransB,
 const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, const double& beta, SDArray<NC>& c)
{
  const size_t K = (NA + NB - NC)/2;
  // calc. contraction shape
  TVector<Dshapes, NC> c_dn_shape;
  gemm_contract_dshape(TransA, TransB, a.dshape(), b.dshape(), c_dn_shape);
  if(c.size() > 0) {
    if(c_dn_shape != c.dshape())
      BTAS_THROW(false, "btas::SDgemm: shape of c mismatched");
    SDscal(beta, c);
  }
  else {
    c.resize(c_dn_shape, false);
  }
  //! call threaded sparse-dgemm
  /*! since c(i, j) = sum_{k} a(i, k) * b(k, j), it's advantageous to store 'b' as b(j, k) order
   *  so that summation over k has less memory-miss when searching non-zero element
   */
  if     (TransA == NoTrans && TransB == NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a, b.transposed_view(K), c);

  else if(TransA == NoTrans && TransB != NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a, b, c);

  else if(TransA != NoTrans && TransB == NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a.transposed_view(K), b.transposed_view(K), c);

  else if(TransA != NoTrans && TransB != NoTrans)
    thread_SDgemm(TransA, TransB, alpha, a.transposed_view(K), b, c);
}

//! Non-BLAS function: (general matrix) * (diagonal matrix)
template<size_t NA, size_t NB>
void SDdimd
(SDArray<NA>& a, const SDArray<NB>& b)
{
  int stride = b.size();
  for(typename SDArray<NA>::iterator ia = a.begin(); ia != a.end(); ++ia) {
    typename SDArray<NB>::const_iterator ib = b.find(ia->first % stride);
    if(ib != b.end()) Ddimd((*ia->second), (*ib->second));
  }
}

//! Non-BLAS function: (diagonal matrix) * (general matrix)
template<size_t NA, size_t NB>
void SDdidm
(const SDArray<NA>& a, SDArray<NB>& b)
{
  const IVector<NB>& b_shape = b.shape();
  int stride = std::accumulate(b_shape.begin()+NA, b_shape.end(), 1, std::multiplies<int>());
  for(typename SDArray<NB>::iterator ib = b.begin(); ib != b.end(); ++ib) {
    typename SDArray<NA>::const_iterator ia = a.find(ib->first / stride);
    if(ia != a.end()) Ddidm((*ia->second), (*ib->second));
  }
}

//! Wrapper function for DBLAS
template<size_t NA, size_t NB, size_t NC>
void SDblasWrapper
(const double& alpha, const SDArray<NA>& a, const SDArray<NB>& b, const double& beta, SDArray<NC>& c)
{
  const size_t K = (NA + NB - NC)/2;
  if(NA == K) {
    // FIXME: this may be wrong or give undesired result
    //        array c needs to be permuted before and after calling Dgemv
    //        otherwise, user has to take care of it
    SDgemv(Trans,   alpha, b, a, beta, c);
  }
  else if(NB == K) {
    SDgemv(NoTrans, alpha, a, b, beta, c);
  }
  else {
    SDgemm(NoTrans, NoTrans, alpha, a, b, beta, c);
  }
}

//! Normalization
template<size_t N>
void SDnormalize(SDArray<N>& x) {
  double norm = SDdot(x, x);
  SDscal(1.0/sqrt(norm), x);
}

//! Orthogonalization
template<size_t N>
void SDorthogonalize(const SDArray<N>& x, SDArray<N>& y) {
  double ovlp = SDdot(x, y);
  SDaxpy(-ovlp, x, y);
}

}; // namespace btas

#endif // _BTAS_CXX11_SDBLAS_H
