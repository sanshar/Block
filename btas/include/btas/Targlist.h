#ifndef _BTAS_ARGMENT_LIST_H
#define _BTAS_ARGMENT_LIST_H 1

#include <vector>
#include <algorithm>
#include <functional>

#ifdef _HAS_INTEL_TBB
#include <tbb/tbb.h>
#endif

#include <btas/btas.h>
#include <btas/DENSE/TArray.h>

namespace btas {

//####################################################################################################
// T_arglist_base
//####################################################################################################

//! Base class of argment list
/*! This class is designed for shared-memory parallelization of sparse-array computations
 *  giving load-balancing functions
 *
 *  TODO: For another design, call() might be implemented as pure-virtual function to carry out
 *  array computation
 */
class T_arglist_base {
protected:
  //! FLOPS count (aprox.) for load-balancing
  size_t m_flops;
public:
  //! Default constructor
  T_arglist_base(size_t flops = 0) : m_flops(flops) { }
  //! Destructor
  virtual ~T_arglist_base() { }
  //! Boolian operator== for compare FLOPS count
  inline bool operator== (const T_arglist_base& other) const { return m_flops == other.m_flops; }
  //! Boolian operator!= for compare FLOPS count
  inline bool operator!= (const T_arglist_base& other) const { return m_flops != other.m_flops; }
  //! Boolian operator< for sorting FLOPS count
  inline bool operator<  (const T_arglist_base& other) const { return m_flops <  other.m_flops; }
  //! Boolian operator> for sorting FLOPS count
  inline bool operator>  (const T_arglist_base& other) const { return m_flops >  other.m_flops; }
};

//####################################################################################################
// R_arglist_* : Argment list for Replication
//####################################################################################################

//! Replication argment list for single argment function
/*! e.g. BLAS SCAL */
template<typename T, size_t N1>
class R_arglist_single : public T_arglist_base {
protected:
  shared_ptr<TArray<T, N1>> m_argment_1;
  //! Function to compute FLOPS count
  virtual size_t mf_flops_count() { return m_argment_1->size(); }
public:
  //! Default constructor
  R_arglist_single() { }
  //! Destructor
  virtual ~R_arglist_single() { }
  //! Initializer
  R_arglist_single
  (const shared_ptr<TArray<T, N1>>& arg1_ptr) {
    m_argment_1 = arg1_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
  //! Reset argment list
  void reset
  (const shared_ptr<TArray<T, N1>>& arg1_ptr) {
    m_argment_1 = arg1_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
};

//! Replication argment list for double argments function
/*! e.g. BLAS COPY, AXPY */
template<typename T, size_t N1, size_t N2>
class R_arglist_double : public T_arglist_base {
protected:
  shared_ptr<TArray<T, N1>> m_argment_1;
  shared_ptr<TArray<T, N2>> m_argment_2;
  //! Function to compute FLOPS count
  virtual size_t mf_flops_count() { return m_argment_1->size(); }
public:
  //! Default constructor
  R_arglist_double() { }
  //! Destructor
  virtual ~R_arglist_double() { }
  //! Initializer
  R_arglist_double
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr) {
    m_argment_1 = arg1_ptr;
    m_argment_2 = arg2_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
  //! Reset argment list
  void reset
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr) {
    m_argment_1 = arg1_ptr;
    m_argment_2 = arg2_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
};

//! Replication argment list for triple argments function
/*! e.g. LAPACK SYEV, HEEV etc. */
template<typename T, size_t N1, size_t N2, size_t N3>
class R_arglist_triple : public T_arglist_base {
protected:
  shared_ptr<TArray<T, N1>> m_argment_1;
  shared_ptr<TArray<T, N2>> m_argment_2;
  shared_ptr<TArray<T, N3>> m_argment_3;
  //! Function to compute FLOPS count
  virtual size_t mf_flops_count() { m_argment_1->size(); }
public:
  //! Default constructor
  R_arglist_triple() { }
  //! Destructor
  virtual ~R_arglist_triple() { }
  //! Initializer
  R_arglist_triple
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr,
   const shared_ptr<TArray<T, N3>>& arg3_ptr) {
    m_argment_1 = arg1_ptr;
    m_argment_2 = arg2_ptr;
    m_argment_3 = arg3_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
  //! Reset argment list
  void reset
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr,
   const shared_ptr<TArray<T, N3>>& arg3_ptr) {
    m_argment_1 = arg1_ptr;
    m_argment_2 = arg2_ptr;
    m_argment_3 = arg3_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
};

//! Replication argment list for quadruple argments function
/*! e.g. LAPACK GESVD etc. */
template<typename T, size_t N1, size_t N2, size_t N3, size_t N4>
class R_arglist_quadra : public T_arglist_base {
protected:
  shared_ptr<TArray<T, N1>> m_argment_1;
  shared_ptr<TArray<T, N2>> m_argment_2;
  shared_ptr<TArray<T, N3>> m_argment_3;
  shared_ptr<TArray<T, N4>> m_argment_4;
  //! Function to compute FLOPS count
  virtual size_t mf_flops_count() { return m_argment_1->size(); }
public:
  //! Default constructor
  R_arglist_quadra() { }
  //! Destructor
  virtual ~R_arglist_quadra() { }
  //! Initializer
  R_arglist_quadra
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr,
   const shared_ptr<TArray<T, N3>>& arg3_ptr,
   const shared_ptr<TArray<T, N4>>& arg4_ptr) {
    m_argment_1 = arg1_ptr;
    m_argment_2 = arg2_ptr;
    m_argment_3 = arg3_ptr;
    m_argment_4 = arg4_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
  //! Reset argment list
  void reset
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr,
   const shared_ptr<TArray<T, N3>>& arg3_ptr,
   const shared_ptr<TArray<T, N4>>& arg4_ptr) {
    m_argment_1 = arg1_ptr;
    m_argment_2 = arg2_ptr;
    m_argment_3 = arg3_ptr;
    m_argment_4 = arg4_ptr;
    T_arglist_base::m_flops = this->mf_flops_count();
  }
};

//####################################################################################################
// C_arglist_* : Argment list for Contruction
//####################################################################################################

//! Argment list for contraction
/*! i.e. BLAS (level 2, 3): GEMV, GER, GEMM, etc.
 *  m_arglist has pairs of array (a, b) and m_c_ptr has the pointer to array (c),
 *  for function call s.t. c = sum_{i} foo(a[i], b[i])
 */
template<typename T, size_t N1, size_t N2, size_t N3>
class C_arglist_triple : public T_arglist_base {
protected:
  std::vector<shared_ptr<TArray<T, N1>>> m_argment_1; //!< array list of a
  std::vector<shared_ptr<TArray<T, N2>>> m_argment_2; //!< array list of b
              shared_ptr<TArray<T, N3>>  m_argment_3; //!< array c
  //! Function to compute FLOPS count
  virtual size_t mf_flops_count
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr) { return m_argment_3 ? m_argment_3->size() : 0; }
public:
  //! Default constructor
  C_arglist_triple() { }
  //! Destructor
  virtual ~C_arglist_triple() { }
  //! Initializer
  C_arglist_triple
  (const shared_ptr<TArray<T, N3>>& arg3_ptr)
  : m_argment_3(arg3_ptr) { }
  //! Reset argment 3
  void reset
  (const shared_ptr<TArray<T, N3>>& arg3_ptr) {
    m_argment_3 = arg3_ptr;
  }
  //! Clear argments
  virtual void clear() {
    m_argment_1.clear();
    m_argment_2.clear();
    m_argment_3.reset();
    T_arglist_base::m_flops = 0;
  }
  //! Add argment list
  void add
  (const shared_ptr<TArray<T, N1>>& arg1_ptr,
   const shared_ptr<TArray<T, N2>>& arg2_ptr) {
    m_argment_1.push_back(arg1_ptr);
    m_argment_2.push_back(arg2_ptr);
    T_arglist_base::m_flops += this->mf_flops_count(arg1_ptr, arg2_ptr);
  }
  //! Return number of argment list
  size_t size() const { return m_argment_1.size(); }
};

//####################################################################################################
// Threaded calling BLAS/LAPACK subroutines
//####################################################################################################

template<class Arglist>
void parallel_call(std::vector<Arglist>& task_list) {
  std::sort(task_list.begin(), task_list.end(), std::greater<Arglist>());
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic) nowait
  for(int i = 0; i < task_list.size(); ++i) {
    task_list[i].call();
  }
}

}; // namespace btas

#endif // _BTAS_ARGMENT_LIST_H
