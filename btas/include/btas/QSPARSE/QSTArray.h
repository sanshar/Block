#ifndef _BTAS_CXX11_QSTARRAY_H
#define _BTAS_CXX11_QSTARRAY_H 1

#include <btas/btas.h>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>

#include <btas/SPARSE/STArray.h>

#include <btas/QSPARSE/Qshapes.h>

namespace btas {

//! Quantum number-based block sparse array
template<typename T, size_t N, class Q = Quantum>
class QSTArray : public STArray<T, N> {
public:
  typedef typename STArray<T, N>::const_iterator const_iterator;
  typedef typename STArray<T, N>::iterator       iterator;

private:
  friend class boost::serialization::access;
  //! Boost serialization
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & boost::serialization::base_object<STArray<T, N>>(*this);
    ar & m_q_total;
    ar & m_q_shape;
  }
  //! Checking non-zero block
  /*! STArray<T, N>::mf_check_allowed is overridden here */
  bool mf_check_allowed(const IVector<N>& block_index) const {
    return (m_q_total == (m_q_shape * block_index));
  }

public:

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Constructors
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Default constructor
  QSTArray() : m_q_total(Q::zero()) { }

  //! Destructor
 ~QSTArray() { }

  //! Construct from quantum number indices
  QSTArray
  (const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape) {
    resize(q_total, q_shape);
  }

  //! Construct from quantum number indices and their dense shapes
  QSTArray(const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape, const TVector<Dshapes, N>& d_shape, bool _allocate = true) {
    resize(q_total, q_shape, d_shape, _allocate);
  }

  //! Construct from quantum number indices and their dense shapes and initialized by constant value
  QSTArray(const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape, const TVector<Dshapes, N>& d_shape, const T& value) {
    resize(q_total, q_shape, d_shape, value);
  }

  //! Construct from quantum number indices and their dense shapes and initialized by gen()
  template<class Generator>
  QSTArray(const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape, const TVector<Dshapes, N>& d_shape, Generator gen) {
    resize(q_total, q_shape, d_shape, gen);
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Copy semantics
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Copy constructor
  QSTArray(const QSTArray& other) : STArray<T, N>(other) {
    m_q_total = other.m_q_total;
    m_q_shape = other.m_q_shape;
  }

  //! Copy assignment operator
  QSTArray& operator= (const QSTArray& other) {
    m_q_total = other.m_q_total;
    m_q_shape = other.m_q_shape;
    STArray<T, N>::copy(other);
    return *this;
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Move semantics
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Move constructor
  QSTArray(QSTArray&& other) : STArray<T, N>(other) {
    m_q_total = std::move(other.m_q_total);
    m_q_shape = std::move(other.m_q_shape);
  }

  //! Move assignment operator
  QSTArray& operator= (QSTArray&& other) {
    m_q_total = std::move(other.m_q_total);
    m_q_shape = std::move(other.m_q_shape);
    STArray<T, N>::operator=(other);
    return *this;
  }

  //! make reference to other
  /*! not complete reference, since elements in m_store are only shared.
   *  so, even if m_shape or m_stride is changed, it won't be affected.
   */
  void reference(const QSTArray& other) {
    m_q_total = other.m_q_total;
    m_q_shape = other.m_q_shape;
    STArray<T, N>::reference(other);
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
  QSTArray subarray(const TVector<Dshapes, N>& _indxs) const {
    QSTArray _ref;
    static_cast<STArray<T, N>&>(_ref) = STArray<T, N>::subarray(_indxs);
    _ref.m_q_total = m_q_total;
    for(int i = 0; i < N; ++i) {
      int nz = _indxs[i].size();
      _ref.m_q_shape[i].resize(nz);
      for(int j = 0; j < nz; ++j)
        _ref.m_q_shape[i][j] = m_q_shape[i].at(_indxs[i][j]);
    }
    return std::move(_ref);
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Resizing functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Construct from quantum number indices
  void resize
  (const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape) {
    IVector<N> s_shape;
    for(int i = 0; i < N; ++i) s_shape[i] = q_shape[i].size();
    m_q_total = q_total;
    m_q_shape = q_shape;
    STArray<T, N>::resize(s_shape);
  }

  //! Resize from quantum number indices and their dense shapes
  void resize
  (const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape, const TVector<Dshapes, N>& d_shape, bool _allocate = true) {
    for(int i = 0; i < N; ++i) assert(q_shape[i].size() == d_shape[i].size());
    m_q_total = q_total;
    m_q_shape = q_shape;
    STArray<T, N>::resize(d_shape, _allocate);
  }

  //! Resize from quantum number indices and their dense shapes and initialized by constant value
  void resize
  (const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape, const TVector<Dshapes, N>& d_shape, const T& value) {
    for(int i = 0; i < N; ++i) assert(q_shape[i].size() == d_shape[i].size());
    m_q_total = q_total;
    m_q_shape = q_shape;
    STArray<T, N>::resize(d_shape, value);
  }

  //! Resize from quantum number indices and their dense shapes and initialized by gen()
  template<class Generator>
  void resize
  (const Q& q_total, const TVector<Qshapes<Q>, N>& q_shape, const TVector<Dshapes, N>& d_shape, Generator gen) {
    for(int i = 0; i < N; ++i) assert(q_shape[i].size() == d_shape[i].size());
    m_q_total = q_total;
    m_q_shape = q_shape;
    STArray<T, N>::resize(d_shape, gen);
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Initializer
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void operator= (const T& value) { STArray<T, N>::operator= (value); }

//template<class Generator>
//void operator= (Generator gen) { STArray<T, N>::operator= (gen); }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Clear and Erase sparse blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Deallocation
  void clear() {
    m_q_total = Q::zero();
    for(int i = 0; i < N; ++i) m_q_shape[i].clear();
    STArray<T, N>::clear();
  }

  //! Erase blocks of which have certain index
  /*! Not good design, since it can't check quantum numbers
   *  \param _rank  rank in which index is associated
   *  \param _index index to be removed
   */
  void erase(int _rank, int _index) {
    STArray<T, N>::erase(_rank, _index);
    m_q_shape[_rank].erase(m_q_shape[_rank].begin()+_index);
  }

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
//  Qshapes<Q> _q_shape;
//  for(size_t i = 0; i < m_shape[_rank]; ++i) {
//    if(_indx_set.find(i) == _indx_set.end())
//      _q_shape.push_back(m_q_shape[_rank][i]);
//  }
//  STArray<T, N>::erase(_rank, _indxs);
//}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Access member variables
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const Q& q() const { return m_q_total; }
  const TVector<Qshapes<Q>, N>& qshape() const { return m_q_shape; }
  const Qshapes<Q>& qshape(int i) const { return m_q_shape[i]; }

  //! Returns quantum numbers corresponding to block-index
  TVector<Q, N> qindex(const IVector<N>& block_index) const {
    return std::move(m_q_shape * block_index);
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parity and Conjugate
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Inplaced parity operation
  /*! For each non-zero block, if q = sum_{i} m_q_shape[p1[i]] has odd particle number, scale -1 */
  void parity(const std::vector<int>& p1) {
    int n = p1.size();
    for(iterator it = this->begin(); it != this->end(); ++it) {
      IVector<N> block_index = this->index(it->first);
      Q qsum = Q::zero();
      for(int i = 0; i < n; ++i)
        qsum = qsum * m_q_shape[p1[i]][block_index[p1[i]]];
      if(qsum.parity())
        it->second->scale(static_cast<T>(-1));
    }
  }

  //! Inplaced parity operation
  /*! For each non-zero block,
   *  scale by sign = prod_{i} ( -1 : both m_q_shape[p1[i]] and m_q_shape[p2[i]] have odd particle number ) */
  void parity(const std::vector<int>& p1, const std::vector<int>& p2) {
    int n = p1.size(); assert(n == p2.size());
    for(iterator it = this->begin(); it != this->end(); ++it) {
      IVector<N> block_index = this->index(it->first);
      bool is_flip_parity = false;
      for(int i = 0; i < n; ++i)
        is_flip_parity ^= m_q_shape[p1[i]][block_index[p1[i]]].parity()
                       && m_q_shape[p2[i]][block_index[p2[i]]].parity();
      if(is_flip_parity)
        it->second->scale(static_cast<T>(-1));
    }
  }

  //! Returns conjugated referece
  /*! Conjugation means flipping direction of quantum indices.  */
  QSTArray conjugate() const {
    QSTArray conj_ref;
    conj_ref.STArray<T, N>::reference(*this);
    conj_ref.m_q_total = -m_q_total;
    for(int i = 0; i < N; ++i)
      conj_ref.m_q_shape[i] = -m_q_shape[i];
    return std::move(conj_ref);
  }

  //! Self conjugation
  void conjugate_self() {
    m_q_total = -m_q_total;
    for(int i = 0; i < N; ++i)
      m_q_shape[i] = -m_q_shape[i];
  }

private:

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Member variables
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //! Total quantum number of array
  Q m_q_total;

  //! Quantum number indices
  /*! m_q_shape[0] = { q[0][0], q[0][1], ... } (for the first  rank)
   *  m_q_shape[1] = { q[1][0], q[1][1], ... } (for the second rank)
   *  ...
   *
   *  e.g. for sparse block A(i, j, k),
   *  if ( m_q_total == m_q_shape[0][i] + m_q_shape[1][j] + m_q_shape[2][k] ) then A(i, j, k) is non-zero */
  TVector<Qshapes<Q>, N>
    m_q_shape;

}; // class QSTArray

}; // namespace btas

template<typename T, size_t N, class Q>
std::ostream& operator<< (std::ostream& ost, const btas::QSTArray<T, N, Q>& a) {
  using std::setw;
  using std::endl;
  ost << "q[T] = " << a.q() << endl;
  for(int i = 0; i < N; ++i)
  ost << "\tq[" << i << "] = " << a.qshape(i) << endl;
  ost << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  ost << static_cast<btas::STArray<T, N>>(a);
  return ost;
}

#endif // _BTAS_CXX11_QSTARRAY_H
