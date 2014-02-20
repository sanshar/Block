#ifndef _BTAS_CXX11_DARGLIST_H
#define _BTAS_CXX11_DARGLIST_H 1

#include <functional>
#include <numeric>

#include <btas/Targlist.h>

#include <btas/DENSE/Dblas.h>
#include <btas/DENSE/Dpermute.h>
#include <btas/DENSE/Dlapack.h>

namespace btas {

//####################################################################################################
// Arglist classes for BLAS (level 1) threaded computation
//####################################################################################################

template<size_t N>
class DcopyArglist : public R_arglist_double<double, N, N> {
public:
  //! Default constructor
  DcopyArglist() { }
  //! Initializer
  DcopyArglist
  (const shared_ptr<DArray<N>>& x_ptr,
   const shared_ptr<DArray<N>>& y_ptr)
  : R_arglist_double<double, N, N>(x_ptr, y_ptr) { }
  //! Call Dcopy
  void call() const { Dcopy(*R_arglist_double<double, N, N>::m_argment_1,
                            *R_arglist_double<double, N, N>::m_argment_2); }
};

template<size_t N>
class DscalArglist : public R_arglist_single<double, N> {
private:
  double
    m_alpha;
public:
  //! Default constructor
  DscalArglist() { }
  //! Initializer
  DscalArglist
  (const double& alpha,
   const shared_ptr<DArray<N>>& x_ptr)
  : R_arglist_single<double, N>(x_ptr), m_alpha(alpha) { }
  //! Reset arglist
  void reset
  (const double& alpha,
   const shared_ptr<DArray<N>>& x_ptr) {
    m_alpha = alpha;
    R_arglist_single<double, N>::reset(x_ptr);
  }
  //! Call Dscal
  void call() const { Dscal(m_alpha, *R_arglist_single<double, N>::m_argment_1); }
};

template<size_t N>
class DaxpyArglist : public R_arglist_double<double, N, N> {
private:
  double
    m_alpha;
public:
  //! Default constructor
  DaxpyArglist() { }
  //! Initializer
  DaxpyArglist
  (const double& alpha,
   const shared_ptr<DArray<N>>& x_ptr,
   const shared_ptr<DArray<N>>& y_ptr)
  : R_arglist_double<double, N, N>(x_ptr, y_ptr), m_alpha(alpha) { }
  //! Reset arglist
  void reset
  (const double& alpha,
   const shared_ptr<DArray<N>>& x_ptr,
   const shared_ptr<DArray<N>>& y_ptr) {
    m_alpha = alpha;
    R_arglist_double<double, N, N>::reset(x_ptr, y_ptr);
  }
  //! Call Daxpy
  void call() const { Daxpy(m_alpha, *R_arglist_double<double, N, N>::m_argment_1,
                                     *R_arglist_double<double, N, N>::m_argment_2); }
};

template<size_t N>
class DpermuteArglist : public R_arglist_double<double, N, N> {
private:
  IVector<N>
    m_permute_index;
public:
  //! Default constructor
  DpermuteArglist() { }
  //! Initializer
  DpermuteArglist
  (const shared_ptr<DArray<N>>& x_ptr,
   const IVector<N>& permute_index,
   const shared_ptr<DArray<N>>& y_ptr)
  : R_arglist_double<double, N, N>(x_ptr, y_ptr),
    m_permute_index(permute_index) { }
  //! Reset arglist
  void reset
  (const shared_ptr<DArray<N>>& x_ptr,
   const IVector<N>& permute_index,
   const shared_ptr<DArray<N>>& y_ptr) {
    m_permute_index = permute_index;
    R_arglist_double<double, N, N>::reset(x_ptr, y_ptr);
  }
  //! Call Dpermute
  void call() const { Dpermute(*R_arglist_double<double, N, N>::m_argment_1, m_permute_index,
                               *R_arglist_double<double, N, N>::m_argment_2); }
};

//####################################################################################################
// Arglist classes for BLAS (level 2, 3) threaded computation
//####################################################################################################

template<size_t NA, size_t NB, size_t NC>
class DgemvArglist : public C_arglist_triple<double, NA, NB, NC> {
private:
  //! this contains index-dependent scaling s.t. parity, Clebsch-Gordan coeff., etc.
  std::vector<double>
    m_scale;
  BTAS_TRANSPOSE
    m_transa;
  double
    m_alpha;
  double
    m_beta;
  //! return FLOPS count
  size_t mf_flops_count
  (const shared_ptr<DArray<NA>>& a_ptr,
   const shared_ptr<DArray<NB>>& b_ptr) { return a_ptr ? a_ptr->size() : 0; }
public:
  //! Default constructor
  DgemvArglist() : m_transa(NoTrans), m_alpha(1.0), m_beta(1.0) { }
  //! Initializer
  DgemvArglist
  (const shared_ptr<DArray<NC>>& c_ptr,
   const BTAS_TRANSPOSE& transa = NoTrans,
   const double& alpha = 1.0,
   const double& beta  = 1.0)
  : C_arglist_triple<double, NA, NB, NC>(c_ptr),
    m_transa(transa), m_alpha(alpha), m_beta(beta) { }
  //! Reset arglist
  void reset
  (const shared_ptr<DArray<NC>>& c_ptr,
   const BTAS_TRANSPOSE& transa = NoTrans,
   const double& alpha = 1.0,
   const double& beta  = 1.0) {
    m_transa = transa;
    m_alpha  = alpha;
    m_beta   = beta;
    C_arglist_triple<double, NA, NB, NC>::reset(c_ptr);
  }
  //! Clear argments
  void clear() {
    m_scale.clear();
    m_transa = NoTrans;
    m_alpha  = 1.0;
    m_beta   = 1.0;
    C_arglist_triple<double, NA, NB, NC>::clear();
  }
  //! Add arglist
  void add
  (const shared_ptr<DArray<NA>>& a_ptr,
   const shared_ptr<DArray<NB>>& b_ptr, double scale = 1.0) {
    m_scale.push_back(scale);
    C_arglist_triple<double, NA, NB, NC>::add(a_ptr, b_ptr);
  }
  //! Call Dgemv
  void call() const {
    for(int i = 0; i < this->size(); ++i)
      Dgemv(m_transa,
            m_scale[i]*m_alpha, *C_arglist_triple<double, NA, NB, NC>::m_argment_1[i],
                                *C_arglist_triple<double, NA, NB, NC>::m_argment_2[i],
                       m_beta,  *C_arglist_triple<double, NA, NB, NC>::m_argment_3);
  }
};

template<size_t NA, size_t NB, size_t NC>
class DgerArglist : public C_arglist_triple<double, NA, NB, NC> {
private:
  //! this contains index-dependent scaling s.t. parity, Clebsch-Gordan coeff., etc.
  std::vector<double>
    m_scale;
  double
    m_alpha;
  //! return FLOPS count
  size_t mf_flops_count
  (const shared_ptr<DArray<NA>>& a_ptr,
   const shared_ptr<DArray<NB>>& b_ptr) { return (a_ptr && b_ptr) ? a_ptr->size()*b_ptr->size() : 0; }
public:
  //! Default constructor
  DgerArglist() : m_alpha(1.0) { }
  DgerArglist
  (const shared_ptr<DArray<NC>>& c_ptr,
   const double& alpha = 1.0)
  : C_arglist_triple<double, NA, NB, NC>(c_ptr),
    m_alpha(alpha) { }
  //! Reset arglist
  void reset
  (const shared_ptr<DArray<NC>>& c_ptr,
   const double& alpha = 1.0) {
    m_alpha = alpha;
    C_arglist_triple<double, NA, NB, NC>::reset(c_ptr);
  }
  //! Clear argments
  void clear() {
    m_scale.clear();
    m_alpha  = 1.0;
    C_arglist_triple<double, NA, NB, NC>::clear();
  }
  //! Add arglist
  void add
  (const shared_ptr<DArray<NA>>& a_ptr,
   const shared_ptr<DArray<NB>>& b_ptr, double scale = 1.0) {
    m_scale.push_back(scale);
    C_arglist_triple<double, NA, NB, NC>::add(a_ptr, b_ptr);
  }
  //! Call Dger
  void call() const {
    for(int i = 0; i < this->size(); ++i)
      Dger(m_scale[i]*m_alpha, *C_arglist_triple<double, NA, NB, NC>::m_argment_1[i],
                               *C_arglist_triple<double, NA, NB, NC>::m_argment_2[i],
                               *C_arglist_triple<double, NA, NB, NC>::m_argment_3);
  }
};

template<size_t NA, size_t NB, size_t NC>
class DgemmArglist : public C_arglist_triple<double, NA, NB, NC> {
private:
  //! this contains index-dependent scaling s.t. parity, Clebsch-Gordan coeff., etc.
  std::vector<double>
    m_scale;
  BTAS_TRANSPOSE
    m_transa;
  BTAS_TRANSPOSE
    m_transb;
  double
    m_alpha;
  double
    m_beta;
  //! return FLOPS count
  size_t mf_flops_count
  (const shared_ptr<DArray<NA>>& a_ptr, const shared_ptr<DArray<NB>>& b_ptr) {
    const size_t K = (NA + NB - NC)/2;
    size_t flops = 0;
    if(a_ptr && b_ptr) {
      flops = a_ptr->size();
      const IVector<NB>& b_shape(b_ptr->shape());
      if(m_transb == NoTrans)
        flops = std::accumulate(b_shape.begin()+K, b_shape.end(), flops, std::multiplies<size_t>());
      else
        flops = std::accumulate(b_shape.begin(), b_shape.begin()+NB-K, flops, std::multiplies<size_t>());
    }
    return flops;
  }
public:
  //! Default constructor
  DgemmArglist() : m_transa(NoTrans), m_transb(NoTrans), m_alpha(1.0), m_beta(1.0) { }
  //! Initializer
  DgemmArglist
  (const shared_ptr<DArray<NC>>& c_ptr,
   const BTAS_TRANSPOSE& transa = NoTrans,
   const BTAS_TRANSPOSE& transb = NoTrans,
   const double& alpha = 1.0,
   const double& beta  = 1.0)
  : C_arglist_triple<double, NA, NB, NC>(c_ptr),
    m_transa(transa), m_transb(transb), m_alpha(alpha), m_beta(beta) { }
  //! Reset arglist
  void reset
  (const shared_ptr<DArray<NC>>& c_ptr,
   const BTAS_TRANSPOSE& transa = NoTrans,
   const BTAS_TRANSPOSE& transb = NoTrans,
   const double& alpha = 1.0,
   const double& beta  = 1.0) {
    m_transa = transa;
    m_transb = transb;
    m_alpha  = alpha;
    m_beta   = beta;
    C_arglist_triple<double, NA, NB, NC>::reset(c_ptr);
  }
  //! Clear argments
  void clear() {
    m_scale.clear();
    m_transa = NoTrans;
    m_transb = NoTrans;
    m_alpha  = 1.0;
    m_beta   = 1.0;
    C_arglist_triple<double, NA, NB, NC>::clear();
  }
  //! Add arglist
  void add
  (const shared_ptr<DArray<NA>>& a_ptr,
   const shared_ptr<DArray<NB>>& b_ptr, double scale = 1.0) {
    m_scale.push_back(scale);
    C_arglist_triple<double, NA, NB, NC>::add(a_ptr, b_ptr);
  }
  //! Call Dgemm
  void call() const {
    for(int i = 0; i < this->size(); ++i)
      Dgemm(m_transa, m_transb,
            m_scale[i]*m_alpha, *C_arglist_triple<double, NA, NB, NC>::m_argment_1[i],
                                *C_arglist_triple<double, NA, NB, NC>::m_argment_2[i],
                       m_beta,  *C_arglist_triple<double, NA, NB, NC>::m_argment_3);
  }
};

//####################################################################################################
// Arglist classes for LAPACK threaded computation
//####################################################################################################

template<size_t NA, size_t NU>
class DgesvdArglist : public R_arglist_quadra<double, NA, 1, NU, NA-NU+2> {
private:
  bool m_calc_u;
  bool m_calc_vt;
public:
  //! Default constructor
  DgesvdArglist() { }
  //! Initializer
  DgesvdArglist
  (const shared_ptr<DArray<NA>     >& a_ptr,
   const shared_ptr<DArray<1>      >& s_ptr,
   const shared_ptr<DArray<NU>     >& u_ptr,
   const shared_ptr<DArray<NA-NU+2>>& v_ptr, bool calc_u = false, bool calc_vt = false)
  : R_arglist_quadra<double, NA, 1, NU, NA-NU+2>(a_ptr, s_ptr, u_ptr, v_ptr), m_calc_u(calc_u), m_calc_vt(calc_vt) { }
  //! Call Dgesvd
  void call() const { Dgesvd(*R_arglist_quadra<double, NA, 1, NU, NA-NU+2>::m_argment_1,
                             *R_arglist_quadra<double, NA, 1, NU, NA-NU+2>::m_argment_2,
                             *R_arglist_quadra<double, NA, 1, NU, NA-NU+2>::m_argment_3,
                             *R_arglist_quadra<double, NA, 1, NU, NA-NU+2>::m_argment_4, m_calc_u, m_calc_vt); }
};

}; // namespace btas

#endif // _BTAS_CXX11_DARGLIST_H
