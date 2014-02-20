#ifndef _BTAS_CXX11_TARRAY_H
#define _BTAS_CXX11_TARRAY_H 1

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include <btas/TVector.h>

namespace btas {

//####################################################################################################
// Forward declaration of TSubArray
//####################################################################################################

template<typename T, size_t N> class TArray;
template<typename T, size_t N> class TSubArray;

//template<typename T, size_t M, size_t N>
//void copy_subarray(const TSubArray<T, M>&, TArray<T, N>&);

//template<typename T, size_t M, size_t N>
//void copy_subarray(const TArray<T, M>&, TSubArray<T, N>&);

/*! \class TArray
 *  \brief Dense array class
 *
 *  Fixed-rank array class implemented in terms of std::array and std::vector
 *  Since using C++11 features, not compatible for C++03 compiler
 *
 *  \param T value type
 *  \param N array rank
 */

template<typename T, size_t N>
class TArray {
private:

  friend class boost::serialization::access;

  //! Any TSubArray classes being friend of TArray<T, N>
  template<typename U, size_t M>
  friend class TSubArray;

  //! Enables to use boost serialization
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_shape & m_stride & m_store; }

public:

  //! TArray<T, N>::iterator
  typedef typename std::vector<T>::iterator       iterator;

  //! TArray<T, N>::const_iterator
  typedef typename std::vector<T>::const_iterator const_iterator;

//####################################################################################################
// Member Functions
//####################################################################################################

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Constructor, Destructor, and Assignment operators
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! default constructor
  TArray() : m_store(new std::vector<T>()) {
    m_shape.fill(0);
    m_stride.fill(0);
  }

  //! destructor
 ~TArray() { }

  //! copy constructor
//explicit TArray(const TArray& other) : m_store(new std::vector<T>()) {
  TArray(const TArray& other) : m_store(new std::vector<T>()) {
    copy(other);
  }

  //! copy assignment operator
  TArray& operator= (const TArray& other) {
    copy(other);
    return *this;
  }

  //! addition assignment operator
  TArray& operator+=(const TArray& other) {
    add (other);
    return *this;
  }

  //! copying from other to this
  void copy(const TArray& other) {
    m_shape  = other.m_shape;
    m_stride = other.m_stride;
    fast_copy(*other.m_store, *m_store);
  }

  //! return copy of this
  TArray copy() const {
    TArray _cpy;
    _cpy.copy(*this);
    return std::move(_cpy);
  }

  //! Copy from sub-array to array
  template<size_t M>
  void copy(const TSubArray<T, M>& a) {
    // Calc. sub-array shape
    IVector<M> a_shape;
    for(int i = 0; i < M; ++i) a_shape[i] = a.m_upper_bound[i]-a.m_lower_bound[i]+1;
    int a_size = std::accumulate(a_shape.begin(), a_shape.end(), 1, std::multiplies<int>());
    assert(m_store->size() == a_size);
    // If 0-dim. array
    if(a_size == 0) return;
    // Striding
    const IVector<M>& a_stride = a.stride();
    int ldt = a_shape[M-1];
    // Get bare pointers
          T* t_ptr = m_store->data();
    const T* a_ptr = a.data();
    // Copying elements
    IVector<M> index(a.m_lower_bound);
    int nrows = a_size / ldt;
    for(int j = 0; j < nrows; ++j, t_ptr += ldt) {
      int offset = dot(a_stride, index);
      _fast_copy(ldt, a_ptr+offset, t_ptr);
      for(int i = static_cast<int>(M)-2; i >= 0; --i) {
        if(++index[i] <= a.m_upper_bound[i]) break;
        index[i] = a.m_lower_bound[i];
      }
    }
  }

  //! scale by const value
  void scale(const T& alpha) {
    fast_scal(alpha, *m_store);
  }

  //! adding  from other to this
  void add (const TArray& other) {
    assert(m_shape  == other.m_shape);
    assert(m_stride == other.m_stride);
    fast_add (*other.m_store, *m_store);
  }

  //! copy constructor from sub-array
  template<size_t M>
  explicit TArray(const TSubArray<T, M>& sub) {
    copy(sub);
  }

  //! copy assignment from sub-array
  template<size_t M>
  TArray& operator= (const TSubArray<T, M>& sub) {
    copy(sub);
    return *this;
  }

  //! move constructor
  explicit TArray(TArray&& other)
  : m_shape(other.m_shape), m_stride(other.m_stride), m_store(other.m_store) { }

  //! move assignment
  TArray& operator= (TArray&& other) {
    m_shape  = std::move(other.m_shape);
    m_stride = std::move(other.m_stride);
    m_store  = std::move(other.m_store);
    return *this;
  }

  //! take data reference from other
  void reference(const TArray& other) {
    m_shape  = other.m_shape;
    m_stride = other.m_stride;
    m_store  = other.m_store;
  }

  //! return data reference of this
  TArray reference() const {
    TArray _ref;
    _ref.reference(*this);
    return std::move(_ref);
  }

  //! convenient constructor with array shape, for N = 1
  explicit TArray(int n01) : m_store(new std::vector<T>()) {
    resize(n01);
  }

  //! convenient constructor with array shape, for N = 2
  TArray(int n01, int n02) : m_store(new std::vector<T>()) {
    resize(n01, n02);
  }

  //! convenient constructor with array shape, for N = 3
  TArray(int n01, int n02, int n03) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03);
  }

  //! convenient constructor with array shape, for N = 4
  TArray(int n01, int n02, int n03, int n04) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04);
  }

  //! convenient constructor with array shape, for N = 5
  TArray(int n01, int n02, int n03, int n04, int n05) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05);
  }

  //! convenient constructor with array shape, for N = 6
  TArray(int n01, int n02, int n03, int n04, int n05, int n06) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05, n06);
  }

  //! convenient constructor with array shape, for N = 7
  TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05, n06, n07);
  }

  //! convenient constructor with array shape, for N = 8
  TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05, n06, n07, n08);
  }

  //! convenient constructor with array shape, for N = 9
  TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05, n06, n07, n08, n09);
  }

  //! convenient constructor with array shape, for N = 10
  TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05, n06, n07, n08, n09, n10);
  }

  //! convenient constructor with array shape, for N = 11
  TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11);
  }

  //! convenient constructor with array shape, for N = 12
  TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11, int n12) : m_store(new std::vector<T>()) {
    resize(n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11, n12);
  }

  //! convenient constructor with array shape, for arbitrary N
  TArray(const IVector<N>& _shape) : m_store(new std::vector<T>()) {
    resize(_shape);
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Resizing functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! resize array shape, for N = 1
  void resize(int n01) {
    IVector< 1> _shape = { n01 };
    resize(_shape);
  }

  //! resize array shape, for N = 2
  void resize(int n01, int n02) {
    IVector< 2> _shape = { n01, n02 };
    resize(_shape);
  }

  //! resize array shape, for N = 3
  void resize(int n01, int n02, int n03) {
    IVector< 3> _shape = { n01, n02, n03 };
    resize(_shape);
  }

  //! resize array shape, for N = 4
  void resize(int n01, int n02, int n03, int n04) {
    IVector< 4> _shape = { n01, n02, n03, n04 };
    resize(_shape);
  }

  //! resize array shape, for N = 5
  void resize(int n01, int n02, int n03, int n04, int n05) {
    IVector< 5> _shape = { n01, n02, n03, n04, n05 };
    resize(_shape);
  }

  //! resize array shape, for N = 6
  void resize(int n01, int n02, int n03, int n04, int n05, int n06) {
    IVector< 6> _shape = { n01, n02, n03, n04, n05, n06 };
    resize(_shape);
  }

  //! resize array shape, for N = 7
  void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07) {
    IVector< 7> _shape = { n01, n02, n03, n04, n05, n06, n07 };
    resize(_shape);
  }

  //! resize array shape, for N = 8
  void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08) {
    IVector< 8> _shape = { n01, n02, n03, n04, n05, n06, n07, n08 };
    resize(_shape);
  }

  //! resize array shape, for N = 9
  void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09) {
    IVector< 9> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09 };
    resize(_shape);
  }

  //! resize array shape, for N = 10
  void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10) {
    IVector<10> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10 };
    resize(_shape);
  }

  //! resize array shape, for N = 11
  void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11) {
    IVector<11> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11 };
    resize(_shape);
  }

  //! resize array shape, for N = 12
  void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11, int n12) {
    IVector<12> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11, n12 };
    resize(_shape);
  }

  //! resize array by _shape, for arbitrary N
  /*! detects rank-mismatching error at compilation time */
  void resize(const IVector<N>& _shape) {
    m_shape = _shape;
    // calculate stride
    size_t stride = 1;
    for(int i = N-1; i >= 0; --i) {
      m_stride[i] = stride;
      stride *= m_shape[i];
    }
    // allocate memory
    m_store->resize(stride);
    return;
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Data Accessing
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! returns first iterator position (const)
  const_iterator begin() const { return m_store->begin(); }

  //! returns first iterator position
        iterator begin()       { return m_store->begin(); }

  //! returns last iterator position (const)
  const_iterator end() const { return m_store->end(); }

  //! returns last iterator position
        iterator end()       { return m_store->end(); }

  //! returns array shape
  const IVector<N>& shape() const { return m_shape; }

  //! returns array shape for rank i
  int shape(int i) const { return m_shape[i]; }

  //! returns array stride
  const IVector<N>& stride() const { return m_stride; }

  //! returns array stride for rank i
  int stride(int i) const { return m_stride[i]; }

  //! returns allocated size
  size_t size() const { return m_store->size(); }

  //! returns array element (N = 1) without range check
  const T& operator() (int i01) const {
    IVector< 1> _index = { i01 };
    return operator()(_index);
  }

  //! returns array element (N = 2) without range check
  const T& operator() (int i01, int i02) const {
    IVector< 2> _index = { i01, i02 };
    return operator()(_index);
  }

  //! returns array element (N = 3) without range check
  const T& operator() (int i01, int i02, int i03) const {
    IVector< 3> _index = { i01, i02, i03 };
    return operator()(_index);
  }

  //! returns array element (N = 4) without range check
  const T& operator() (int i01, int i02, int i03, int i04) const {
    IVector< 4> _index = { i01, i02, i03, i04 };
    return operator()(_index);
  }

  //! returns array element (N = 5) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05) const {
    IVector< 5> _index = { i01, i02, i03, i04, i05 };
    return operator()(_index);
  }

  //! returns array element (N = 6) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06) const {
    IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
    return operator()(_index);
  }

  //! returns array element (N = 7) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07) const {
    IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
    return operator()(_index);
  }

  //! returns array element (N = 8) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) const {
    IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
    return operator()(_index);
  }

  //! returns array element (N = 9) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) const {
    IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
    return operator()(_index);
  }

  //! returns array element (N = 10) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) const {
    IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
    return operator()(_index);
  }

  //! returns array element (N = 11) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) const {
    IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
    return operator()(_index);
  }

  //! returns array element (N = 12) without range check
  const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) const {
    IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
    return operator()(_index);
  }

  //! returns array element (arbitrary N) without range check
  const T& operator() (const IVector<N>& _index) const {
    return (*m_store)[dot(_index, m_stride)];
  }

  //! returns array element (N = 1) without range check
  T& operator() (int i01) {
    IVector< 1> _index = { i01 };
    return operator()(_index);
  }

  //! returns array element (N = 2) without range check
  T& operator() (int i01, int i02) {
    IVector< 2> _index = { i01, i02 };
    return operator()(_index);
  }

  //! returns array element (N = 3) without range check
  T& operator() (int i01, int i02, int i03) {
    IVector< 3> _index = { i01, i02, i03 };
    return operator()(_index);
  }

  //! returns array element (N = 4) without range check
  T& operator() (int i01, int i02, int i03, int i04) {
    IVector< 4> _index = { i01, i02, i03, i04 };
    return operator()(_index);
  }

  //! returns array element (N = 5) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05) {
    IVector< 5> _index = { i01, i02, i03, i04, i05 };
    return operator()(_index);
  }

  //! returns array element (N = 6) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05, int i06) {
    IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
    return operator()(_index);
  }

  //! returns array element (N = 7) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07) {
    IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
    return operator()(_index);
  }

  //! returns array element (N = 8) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) {
    IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
    return operator()(_index);
  }

  //! returns array element (N = 9) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) {
    IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
    return operator()(_index);
  }

  //! returns array element (N = 10) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) {
    IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
    return operator()(_index);
  }

  //! returns array element (N = 11) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) {
    IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
    return operator()(_index);
  }

  //! returns array element (N = 12) without range check
  T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) {
    IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
    return operator()(_index);
  }

  //! returns array element (arbitrary N) without range check
  T& operator() (const IVector<N>& _index) {
    return (*m_store)[dot(_index, m_stride)];
  }

  //! returns array element (N = 1) with range check
  const T& at(int i01) const {
    IVector< 1> _index = { i01 };
    return at(_index);
  }

  //! returns array element (N = 2) with range check
  const T& at(int i01, int i02) const {
    IVector< 2> _index = { i01, i02 };
    return at(_index);
  }

  //! returns array element (N = 3) with range check
  const T& at(int i01, int i02, int i03) const {
    IVector< 3> _index = { i01, i02, i03 };
    return at(_index);
  }

  //! returns array element (N = 4) with range check
  const T& at(int i01, int i02, int i03, int i04) const {
    IVector< 4> _index = { i01, i02, i03, i04 };
    return at(_index);
  }

  //! returns array element (N = 5) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05) const {
    IVector< 5> _index = { i01, i02, i03, i04, i05 };
    return at(_index);
  }

  //! returns array element (N = 6) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05, int i06) const {
    IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
    return at(_index);
  }

  //! returns array element (N = 7) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07) const {
    IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
    return at(_index);
  }

  //! returns array element (N = 8) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) const {
    IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
    return at(_index);
  }

  //! returns array element (N = 9) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) const {
    IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
    return at(_index);
  }

  //! returns array element (N = 10) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) const {
    IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
    return at(_index);
  }

  //! returns array element (N = 11) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) const {
    IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
    return at(_index);
  }

  //! returns array element (N = 12) with range check
  const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) const {
    IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
    return at(_index);
  }

  //! returns array element (arbitrary N) with range check
  const T& at(const IVector<N>& _index) const {
    return m_store->at(dot(_index, m_stride));
  }

  //! returns array element (N = 1) with range check
  T& at(int i01) {
    IVector< 1> _index = { i01 };
    return at(_index);
  }

  //! returns array element (N = 2) with range check
  T& at(int i01, int i02) {
    IVector< 2> _index = { i01, i02 };
    return at(_index);
  }

  //! returns array element (N = 3) with range check
  T& at(int i01, int i02, int i03) {
    IVector< 3> _index = { i01, i02, i03 };
    return at(_index);
  }

  //! returns array element (N = 4) with range check
  T& at(int i01, int i02, int i03, int i04) {
    IVector< 4> _index = { i01, i02, i03, i04 };
    return at(_index);
  }

  //! returns array element (N = 5) with range check
  T& at(int i01, int i02, int i03, int i04, int i05) {
    IVector< 5> _index = { i01, i02, i03, i04, i05 };
    return at(_index);
  }

  //! returns array element (N = 6) with range check
  T& at(int i01, int i02, int i03, int i04, int i05, int i06) {
    IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
    return at(_index);
  }

  //! returns array element (N = 7) with range check
  T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07) {
    IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
    return at(_index);
  }

  //! returns array element (N = 8) with range check
  T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) {
    IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
    return at(_index);
  }

  //! returns array element (N = 9) with range check
  T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) {
    IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
    return at(_index);
  }

  //! returns array element (N = 10) with range check
  T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) {
    IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
    return at(_index);
  }

  //! returns array element (N = 11) with range check
  T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) {
    IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
    return at(_index);
  }

  //! returns array element (N = 12) with range check
  T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) {
    IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
    return at(_index);
  }

  //! returns array element (arbitrary N) with range check
  T& at(const IVector<N>& _index) {
    return m_store->at(dot(_index, m_stride));
  }

  //! slice array to return sub-array object
  /*! sub-array is constructed from array elements [lbound[0]:ubound[0], lbound[1]:ubound[1], ...] */
  TSubArray<T, N> subarray(const IVector<N>& lbound, const IVector<N>& ubound) const {
    return TSubArray<T, N>(*this, lbound, ubound);
  }

  //! returns the first pointer of array elements
  const T* data() const { return m_store->data(); }

  //! returns the first pointer of array elements
        T* data()       { return m_store->data(); }

  //! fills elements by constant value
  void fill(const T& val) {
    std::fill(m_store->begin(), m_store->end(), val);
  }

  //! fills elements by constant value
  void operator= (const T& val) { fill(val); }

  //! generates array elements by function gen
  /*! Generator is either default constructible class or function pointer, which can be called by gen() */
  template<class Generator>
  void generate(Generator gen) {
    std::generate(m_store->begin(), m_store->end(), gen);
  }

////! generates array elements by function gen
//template<class Generator>
//void operator= (Generator gen) { generate(gen); }

  //! deallocate storage
  void clear() {
    m_shape.fill(0);
    m_stride.fill(0);
    m_store->clear();
  }

private:

//####################################################################################################
// Member Variables
//####################################################################################################

  //! array shape
  IVector<N>
    m_shape;

  //! array stride
  IVector<N>
    m_stride;

  //! array storage
  shared_ptr<std::vector<T>>
    m_store;

}; // class TArray

}; // namespace btas

#include <btas/DENSE/TSubArray.h>

#endif // _BTAS_CXX11_TARRAY_H
