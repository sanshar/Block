#ifndef _BTAS_CXX11_STARRAY_H
#define _BTAS_CXX11_STARRAY_H 1

#include <iostream>
#include <iomanip>
#include <map>
#include <memory>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <btas/btas.h>
#include <btas/TVector.h>

#include <btas/DENSE/TArray.h>

namespace btas {

//! Block-sparse array class
/*!
 *  explain here for more detail
 */

template<typename T, size_t N>
class STArray {
private:
  // Alias to data type
  typedef std::map<int, shared_ptr<TArray<T, N>>> DataType;

public:
  // Alias to iterator
  typedef typename DataType::const_iterator const_iterator;
  typedef typename DataType::iterator       iterator;

private:
  // Boost serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_shape & m_dn_shape & m_stride & m_store; }

protected:
  //! Checking non-zero block
  /*! This should be overridden so that non-zero block can be determined from STArray class */
  virtual bool mf_check_allowed(const IVector<N>& _index) const { return true; }

  // KEEP FOR INSERTION CHECK

  void mf_check_dshape(const IVector<N>& _index, const IVector<N>& _shape) {
    IVector<N> _chk_shape = m_dn_shape & _index;
    for(int i = 0; i < N; ++i) {
      if(_chk_shape[i] == _shape[i]) continue;
      BTAS_THROW(_chk_shape[i] == 0, "btas::STArray::mf_check_dshape: requested shape is inconsistent");
      m_dn_shape[i][_index[i]] = _shape[i]; // update dense-block shape
    }
  }

public:

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Constructors
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Default constructor
  STArray() { m_shape.fill(0); m_stride.fill(0); }

  //! Destructor
  virtual ~STArray() { }

  //! Construct by sparse-block shape
  STArray(const IVector<N>& _shape) { resize(_shape); }

  //! Construct by dense-block shapes
  STArray(const TVector<Dshapes, N>& _dn_shape, bool _allocate = true) { resize(_dn_shape, _allocate); }

  //! Construct by dense-block shapes and fill elements by value
  STArray(const TVector<Dshapes, N>& _dn_shape, const T& value) { resize(_dn_shape, value); }

  //! Construct by dense-block shapes and fill elements by gen()
  template<class Generator>
  STArray(const TVector<Dshapes, N>& _dn_shape, Generator gen) { resize(_dn_shape, gen); }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy semantics
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Copy constructor
  STArray(const STArray& other) { copy(other); }

  //! Copy assignment operator
  STArray& operator= (const STArray& other) { copy(other); return *this; }

  //! Take deep copy of other
  void copy(const STArray& other) {
    m_shape    = other.m_shape;
    m_dn_shape = other.m_dn_shape;
    m_stride   = other.m_stride;
    m_store.clear();
    iterator ip = m_store.begin();
    for(const_iterator it = other.m_store.begin(); it != other.m_store.end(); ++it) {
      if(!it->second) continue; // remove NULL element upon copying
      ip = m_store.insert(ip, std::make_pair(it->first, shared_ptr<TArray<T, N>>(new TArray<T, N>(*(it->second)))));
    }
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Move semantics
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Move constructor
  STArray(STArray&& other) {
    m_shape     = std::move(other.m_shape);
    m_dn_shape  = std::move(other.m_dn_shape);
    m_stride    = std::move(other.m_stride);
    m_store     = std::move(other.m_store);
  }

  //! Move assignment operator
  STArray& operator= (STArray&& other) {
    m_shape     = std::move(other.m_shape);
    m_dn_shape  = std::move(other.m_dn_shape);
    m_stride    = std::move(other.m_stride);
    m_store     = std::move(other.m_store);
    return *this;
  }

  //! make reference to other
  /*! not complete reference, since elements in m_store are only shared.
   *  so, even if m_shape or m_stride is changed, it won't be affected.
   */
  void reference(const STArray& other) {
    m_shape     = other.m_shape;
    m_dn_shape  = other.m_dn_shape;
    m_stride    = other.m_stride;
    m_store     = other.m_store;
  }

  //! Make subarray reference
  /*! \param _indxs contains subarray indices
   *  e.g.
   *  sparse shape = { 4, 4 }
   *  _indxs = { { 1, 3 }, { 0, 2, 3} }
   *
   *     0  1  2  3           0  2  3
   *    +--+--+--+--+        +--+--+--+
   *  0 |  |  |  |  |  ->  1 |**|**|**|
   *    +--+--+--+--+        +--+--+--+
   *  1 |**|  |**|**|      3 |**|**|**|
   *    +--+--+--+--+        +--+--+--+
   *  2 |  |  |  |  |
   *    +--+--+--+--+
   *  3 |**|  |**|**|
   *    +--+--+--+--+
   *
   *  ** blocks are only kept to make subarray
   */
  STArray subarray(const TVector<Dshapes, N>& _indxs) const {
    TVector<Dshapes, N> _indx_map;
    TVector<Dshapes, N> _dn_shape;
    for(int i = 0; i < N; ++i) {
      _indx_map[i].resize(m_shape[i], -1);
      _dn_shape[i].reserve(_indxs[i].size());
      int n = 0;
      for(int j = 0; j < _indxs[i].size(); ++j) {
        _indx_map[i].at(_indxs[i][j]) = n++;
        _dn_shape[i].push_back(m_dn_shape[i].at(_indxs[i][j]));
      }
    }

    STArray _ref(_dn_shape, false);

    iterator ip = _ref.m_store.begin();
    for(const_iterator it = m_store.begin(); it != m_store.end(); ++it) {
      IVector<N> _index = _indx_map & index(it->first);
      bool kept = true;
      for(int i = 0; kept && i < N; ++i) {
        kept &= (_index[i] >= 0);
      }
      if(kept) ip = _ref.m_store.insert(ip, std::make_pair(_ref.tag(_index), it->second));
    }
    return std::move(_ref);
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Resizing functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Resize by sparse-block shape
  void resize(const IVector<N>& _shape) {
    m_shape = _shape;
    for(int i = 0; i < N; ++i) {
      m_dn_shape[i] = Dshapes(m_shape[i], 0);
    }
    int stride = 1;
    for(int i = N-1; i >= 0; --i) {
      m_stride[i] = stride;
      stride *= m_shape[i];
    }
    m_store.clear();
  }

  //! Resize by dense-block shapes using this->mf_check_allowed(index)
  void resize(const TVector<Dshapes, N>& _dn_shape, bool _allocate = true) {
    // calc. sparse-block shape
    IVector<N> _shape;
    for(int i = 0; i < N; ++i) _shape[i] = _dn_shape[i].size();
    resize(_shape);
    m_dn_shape = _dn_shape;

    if(_allocate) allocate();
  }

  //! Allocate all allowed blocks (existed blocks are collapsed)
  void allocate() {
    iterator it = m_store.begin();
    IVector<N> _index = uniform<int, N>(0);

    m_store.clear();

    for(size_t ib = 0; ib < size(); ++ib) {
      // assume derived mf_check_allowed being called
      if(this->mf_check_allowed(_index)) {
        if(m_dn_shape * _index > 0) // check non-zero size
          it = m_store.insert(it, std::make_pair(ib, shared_ptr<TArray<T, N>>(new TArray<T, N>(m_dn_shape & _index))));
      }
      // index increment
      for(int id = N-1; id >= 0; --id) {
        if(++_index[id] < m_shape[id]) break;
        _index[id] = 0;
      }
    }
  }

  //! Resize by dense-block shapes and fill all elements by value
  void resize(const TVector<Dshapes, N>& _dn_shape, const T& value) { resize(_dn_shape); fill(value); }

  //! Resize by dense-block shapes and fill all elements by gen()
  template<class Generator>
  void resize(const TVector<Dshapes, N>& _dn_shape, Generator gen) { resize(_dn_shape); generate(gen); }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Fill and Generage elements
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! fills all elements by value
  void fill(const T& value) {
    for(iterator it = m_store.begin(); it != m_store.end(); ++it) it->second->fill(value);
  }

  //! fills all elements by value
  void operator= (const T& value) { fill(value); }

  //! generates all elements by gen()
  template<class Generator>
  void generate(Generator gen) {
    for(iterator it = m_store.begin(); it != m_store.end(); ++it) it->second->generate(gen);
  }

////! generates all elements by gen()
//template<class Generator>
//void operator= (Generator gen) { generate(gen); } // this is ambiguous with copy assignment operaotr

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Clear and Erase sparse blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Erase certain block by index
  void erase (const IVector<N>& _index) { m_store.erase(tag(_index)); }

  //! Erase certain block by tag
  void erase (const int& _tag) {  m_store.erase(_tag); }

  //! Deallocation
  virtual void clear() {
    m_shape.fill(0);
    m_stride.fill(0);
    for(int i = 0; i < N; ++i) m_dn_shape[i].clear();
    m_store.clear();
  }

  //! Erase blocks of which have certain index
  /*!
   *  \param _rank  rank in which index is associated
   *  \param _index index to be removed
   */
  virtual void erase(int _rank, int _index_erase) {
    assert(_index_erase >= 0 && _index_erase < m_shape[_rank]);
    IVector<N> _shape = m_shape; --_shape[_rank];
    STArray _ref(_shape);
    iterator ip = _ref.m_store.begin();
    for(iterator it = m_store.begin(); it != m_store.end(); ++it) {
      IVector<N> _index = index(it->first);
      if(_index[_rank] == _index_erase) continue;
      if(_index[_rank] >  _index_erase) --_index[_rank];
      ip = _ref.m_store.insert(ip, std::make_pair(_ref.tag(_index), it->second));
    }
    *this = std::move(_ref);
  } //				FIXME: This might not be good. Duplicated with subarray function

////! Erase blocks of which have certain set of indices
///*!
// *  \param _rank  rank in which index is associated
// *  \param _indxs indices to be removed
// */
//virtual void erase(int _rank, const std::vector<int>& _indxs) {
//  std::set<int> _indx_set(_indxs.begin(), _indxs.end());
//  assert(_indx_set.size() <= m_shape[_rank]);
//  assert(*_indx_set.begin() >= 0 && *_indx_set.rbegin() < m_shape[_rank]);
//  std::map<int, int> _indx_map;
//  int nnz = 0;
//  for(size_t i = 0; i < m_shape[_rank]; ++i) {
//    if(_indx_set.find(i) == _indx_set.end())
//      _indx_map.insert(std::make_pair(i, nnz++));
//  }
//  IVector<N> _shape = m_shape; _shape[_rank] = nnz;
//  STArray _ref(_shape);
//  iterator ip = _ref.m_store.begin();
//  for(iterator it = m_store.begin(); it != m_store.end(); ++it) {
//    IVector<N> _index = indxs(it->first);
//    auto imap = _indx_map.find(_index[_rank]);
//    if(imap == _indx_map.end()) continue;
//    _index[_rank] = imap->second;
//    ip = _ref.m_store.insert(ip, std::make_pair(_ref.tag(_index), it->second));
//  }
//  *this = std::move(_ref);
//}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Index <--> Tag conversion
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! convert tag to index
  IVector<N> index(int _tag) const {
    IVector<N> _index;
    for(int i = 0; i < N; ++i) {
      _index[i] = _tag / m_stride[i];
      _tag      = _tag % m_stride[i];
    }
    return std::move(_index);
  }

  //! convert index to tag
  int tag(const IVector<N>& _index) const { return dot(_index, m_stride); }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Access member variables
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Returns sparse-block shape
  const IVector<N>& shape() const { return m_shape; }

  //! Returns sparse-block shape for rank i
  const int& shape(int i) const { return m_shape[i]; }

  //! Returns sparse-block stride
  const IVector<N>& stride() const { return m_stride; }

  //! Returns sparse-block stride for rank i
  const int& stride(int i) const { return m_stride[i]; }

  //! Returns number of non-zero sparse-blocks
  size_t nnz() const { return m_store.size(); }

  //! Returns total number of sparse-blocks (includes zero blocks)
  size_t size() const { return m_stride[0]*m_shape[0]; }

  //! Returns dense-block shapes
  const TVector<Dshapes, N>& dshape() const { return m_dn_shape; }

  //! Returns dense-block shapes for rank i
  const Dshapes& dshape(int i) const { return m_dn_shape[i]; }

  //! Check and update dense-block shapes
  const TVector<Dshapes, N>& check_dshape() {
    for(const_iterator it = m_store.begin(); it != m_store.end(); ++it)
      mf_check_dshape(index(it->first), it->second->shape());
    return m_dn_shape;
  }

  //! Calc. and return net dense-block shapes
  TVector<Dshapes, N> check_net_dshape() const {
    TVector<Dshapes, N> _net_dshape;
    for(size_t i = 0; i < N; ++i) _net_dshape[i].resize(m_shape[i], 0);
    for(const_iterator it = m_store.begin(); it != m_store.end(); ++it) {
      IVector<N> _index = index(it->first);
      for(int i = 0; i < N; ++i) {
        if(_net_dshape[i][_index[i]] > 0) {
          BTAS_THROW(_net_dshape[i][_index[i]] == it->second->shape(i), "btas::STArray::check_net_dshape: found mismatched dense-block shape");
        }
        else {
          _net_dshape[i][_index[i]] = it->second->shape(i);
        }
      }
    }
    return std::move(_net_dshape);
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Iterators: Definitions are related to std::map
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const_iterator begin() const { return m_store.begin(); }
        iterator begin()       { return m_store.begin(); }

  const_iterator end() const { return m_store.end(); }
        iterator end()       { return m_store.end(); }

  const_iterator find(const IVector<N>& _index) const { return m_store.find(tag(_index)); }
        iterator find(const IVector<N>& _index)       { return m_store.find(tag(_index)); }

  const_iterator lower_bound(const IVector<N>& _index) const { return m_store.lower_bound(tag(_index)); }
        iterator lower_bound(const IVector<N>& _index)       { return m_store.lower_bound(tag(_index)); }

  const_iterator upper_bound(const IVector<N>& _index) const { return m_store.upper_bound(tag(_index)); }
        iterator upper_bound(const IVector<N>& _index)       { return m_store.upper_bound(tag(_index)); }

  const_iterator find(const int& _tag) const { return m_store.find(_tag); }
        iterator find(const int& _tag)       { return m_store.find(_tag); }

  const_iterator lower_bound(const int& _tag) const { return m_store.lower_bound(_tag); }
        iterator lower_bound(const int& _tag)       { return m_store.lower_bound(_tag); }

  const_iterator upper_bound(const int& _tag) const { return m_store.upper_bound(_tag); }
        iterator upper_bound(const int& _tag)       { return m_store.upper_bound(_tag); }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Insert dense-block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! return true if the requested block is non-zero, called by block tag
  bool allowed(const int& _tag) const { return this->mf_check_allowed(index(_tag)); }
  //! return true if the requested block is non-zero, called by block index
  bool allowed(const IVector<N>& _index) const { return this->mf_check_allowed(_index); }

  //! reserve non-zero block and return its iterator, by block tag
  /*! if the requested block already exists:
   *  - return its iterator
   *  - or, return error if it's not allowed
   *  if the requested block hasn't allocated
   *  - allocate dense-array block and return its iterator
   *  - or, return last iterator if it's not allowed, with warning message (optional)
   */
  iterator reserve(const int& _tag) {
    IVector<N> _index = index(_tag);
    // check if the requested block can be non-zero
    iterator it = find(_tag);
    if(this->mf_check_allowed(_index)) {
      if(it == end()) {
        it = m_store.insert(it, std::make_pair(_tag, shared_ptr<TArray<T, N>>(new TArray<T, N>(m_dn_shape & _index))));
      }
      else {
        BTAS_THROW((m_dn_shape & _index) == it->second->shape(), "btas::STArray::reserve: existed block has inconsistent shape");
      }
    }
    else {
      if(it != end()) {
        BTAS_THROW(false, "btas::STArray::reserve: non-zero block already exists despite it must be zero");
      }
#ifdef _PRINT_WARNINGS
      else {
        BTAS_DEBUG("WARNING: btas::STArray::reserve: requested block must be zero, returns end()");
      }
#endif
    }
    return it;
  }

  //! reserve non-zero block and return its iterator, by block index
  iterator reserve(const IVector<N>& _index) {
    int _tag = tag(_index);
    // check if the requested block can be non-zero
    iterator it = find(_tag);
    if(this->mf_check_allowed(_index)) {
      if(it == end()) {
        it = m_store.insert(it, std::make_pair(_tag, shared_ptr<TArray<T, N>>(new TArray<T, N>(m_dn_shape & _index))));
      }
      else {
        BTAS_THROW((m_dn_shape & _index) == it->second->shape(), "btas::STArray::reserve: existed block has inconsistent shape");
      }
    }
    else {
      if(it != end()) {
        BTAS_THROW(false, "btas::STArray::reserve; non-zero block already exists despite it must be zero");
      }
#ifdef _PRINT_WARNINGS
      else {
        BTAS_DEBUG("WARNING: btas::STArray::reserve: requested block must be zero, returns end()");
      }
#endif
    }
    return it;
  }

  //! insert dense-array block and return its iterator, by block tag
  /*! if the requested block already exists:
   *  - add array to it, return its iterator
   *  if the requested block hasn't allocated
   *  - insert dense-array block and return its iterator
   *  - or, return last iterator if it's not allowed, with warning message (optional)
   */
  iterator insert(const int& _tag, const TArray<T, N>& block) {
    IVector<N> _index = index(_tag);
    // check if the requested block can be non-zero
    iterator it = m_store.end();
    if(this->mf_check_allowed(_index)) {
      mf_check_dshape(_index, block.shape());
      it = find(_tag);
      if(it != end())
        it->second->add(block);
      else
        it = m_store.insert(it, std::make_pair(_tag, shared_ptr<TArray<T, N>>(new TArray<T, N>(block))));
    }
#ifdef _PRINT_WARNINGS
    else {
      BTAS_DEBUG("WARNING: btas::STArray::insert: requested block must be zero, unable to be inserted");
    }
#endif
    return it;
  }

  //! insert dense-array block and return its iterator, by block index
  iterator insert(const IVector<N>& _index, const TArray<T, N>& block) {
    int _tag = tag(_index);
    // check if the requested block can be non-zero
    iterator it = m_store.end();
    if(this->mf_check_allowed(_index)) {
      mf_check_dshape(_index, block.shape());
      it = find(_tag);
      if(it != end())
        it->second->add(block);
      else
        it = m_store.insert(it, std::make_pair(_tag, shared_ptr<TArray<T, N>>(new TArray<T, N>(block))));
    }
#ifdef _PRINT_WARNINGS
    else {
      BTAS_DEBUG("WARNING: btas::STArray::insert: requested block must be zero, unable to be inserted");
    }
#endif
    return it;
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Transposed and Permuted references
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! return reference in which the sparse-blocks are transposed
  /*! dense-arrays are not transposed
   *  \param K is a rank to be transposed
   *  e.g. N = 6, K = 4
   *  [i,j,k,l,m,n] -> [m,n,i,j,k,l]
   *  |1 2 3 4|5 6|    |5 6|1 2 3 4|
   */
  STArray transposed_view(int K) const {
    STArray trans;
    if(K == N) {
      trans.reference(*this);
    }
    else {
      TVector<Dshapes, N> t_dn_shape = transpose(m_dn_shape, K);
      trans.resize(t_dn_shape, false);

      int oldstr = m_stride[K-1];
      int newstr = size() / oldstr;
      iterator ip = trans.m_store.begin();
      for(const_iterator it = m_store.begin(); it != m_store.end(); ++it) {
        int oldtag = it->first;
        int newtag = oldtag / oldstr + (oldtag % oldstr)*newstr;
        ip = trans.m_store.insert(ip, std::make_pair(newtag, it->second));
      }
    }
    return std::move(trans);
  }

  //! return reference in which the sparse-blocks are permuted by pindex
  /*! dense-arrays are not permuted */
  STArray permuted_view(const IVector<N>& pindex) const {
    STArray pmute;
    if(pindex == sequence<N>(0, 1)) {
      pmute.reference(*this);
    }
    else {
      TVector<Dshapes, N> p_dn_shape = permute(m_dn_shape, pindex);
      pmute.resize(p_dn_shape, false);

      IVector<N> p_stride;
      for(int i = 0; i < N; ++i) p_stride[pindex[i]] = pmute.m_stride[i];
      iterator ip = pmute.m_store.begin();
      for(const_iterator it = m_store.begin(); it != m_store.end(); ++it) {
        IVector<N> _index(index(it->first));
        ip = pmute.m_store.insert(ip, std::make_pair(dot(_index, p_stride), it->second));
      }
    }
    return pmute;
  }

protected:

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Member variables
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! sparse-block shape
  IVector<N>
    m_shape;

  //! dense-block shapes
  TVector<Dshapes, N>
    m_dn_shape;

  //! stride for sparse-block
  IVector<N>
    m_stride;

  //! non-zero data array mapped by tag
  DataType
    m_store;

}; // class STArray

}; // namespace btas

//! C++ style printing function
template<typename T, size_t N>
std::ostream& operator<< (std::ostream& ost, const btas::STArray<T, N>& a) {
  using std::endl;
  // print out sparsity information
  const btas::IVector<N>& a_shape = a.shape();
  ost << "block shape = [ ";
  for(int i = 0; i < N-1; ++i) ost << a_shape[i] << " x ";
  ost << a_shape[N-1] << " ] ( sparsity = " << a.nnz() << " / " << a.size() << " ) " << endl;

  for(typename btas::STArray<T, N>::const_iterator ib = a.begin(); ib != a.end(); ++ib) {
    ost << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    btas::IVector<N> b_index = a.index(ib->first);
    ost << "\tindex = [ ";
    for(int i = 0; i < N-1; ++i) ost << b_index[i] << ", ";
    ost << b_index[N-1] << " ] : " << *ib->second << endl;
  }
  return ost;
}

#endif // _BTAS_CXX11_STARRAY_H
