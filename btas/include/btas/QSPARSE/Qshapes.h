#ifndef _BTAS_CXX11_QSHAPES_H
#define _BTAS_CXX11_QSHAPES_H 1

#include <vector>
#include <set>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

// If user doesn't have default quantum class, BTAS provide the default class
#ifdef _DEFAULT_QUANTUM
#include <btas/QSPARSE/Quantum.h>
#endif

#include <btas/TVector.h>

namespace btas {

//! Quantum number vector
/*! Since typedef of std::vector<Q> is ambiguous to define operators,
 *  it's implemented as a different class in terms of std::vector<Q>.
 *
 *  class Q must have overloaded operators: =, + (unary), - (unary), *, ==, !=, <, >
 *  It's also required,
 *  - const static function Q::zero() which gives zero quantum number,
 *  - boolian function parity() gives true when it has odd particle number in the case of fermion
 *  - and, clebsch() which gives Clebsch-Gordan coefficient in the case non-Abelian symmetry
 *
 *  Default quantum number class is named by 'btas::Quantum'
 */
template<class Q = Quantum>
class Qshapes : public std::vector<Q> {
private:
  friend class boost::serialization::access;
  //! Boost serialization
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & boost::serialization::base_object<std::vector<Q>>(*this);
  }
public:
  typedef typename std::vector<Q>::value_type value_type;
  typedef typename std::vector<Q>::size_type size_type;
  typedef typename std::vector<Q>::iterator iterator;
  typedef typename std::vector<Q>::const_iterator const_iterator;
  //  Constructors in std::vector
  Qshapes() { }
  explicit Qshapes(size_type n)                    : std::vector<Q>(n)           { }
  Qshapes(size_type n, const value_type& val)      : std::vector<Q>(n, val)      { }
  template<class InputIterator>
  Qshapes(InputIterator first, InputIterator last) : std::vector<Q>(first, last) { }
  Qshapes(const Qshapes& x)                        : std::vector<Q>(x)           { }
  Qshapes(Qshapes&& x)                             : std::vector<Q>(x)           { }
  Qshapes(std::initializer_list<value_type> il)    : std::vector<Q>(il)          { }
 ~Qshapes() { }
  //  Overloaded operators
  //! Copy assignment operator
  Qshapes& operator= (const Qshapes&  other) { std::vector<Q>::operator= (other); return *this; }
  //! Move assignment operator
  Qshapes& operator= (      Qshapes&& other) { std::vector<Q>::operator= (other); return *this; }
  //! Contract two quantum numbers: { q(ij) : q(i) * q(j) }
  Qshapes  operator* (const Qshapes& other) const {
    Qshapes qij;
    qij.reserve(this->size()*other.size());
    for(const Q& qi : *this)
      for(const Q& qj : other) qij.push_back(qi * qj);
    return std::move(qij);
  }
  //! Contract two quantum numbers and chose unique quantum numbers: { q(ij) : q(i) * q(j) } => { q(st) }
  /*! e.g. { q(ij) = q1, q2, q1, q3 } => { q(st) = q1, q2, q3 } */
  Qshapes  operator& (const Qshapes& other) const {
    std::set<Q> qij;
    for(const Q& qi : *this)
      for(const Q& qj : other) qij.insert(qi * qj);
    Qshapes qst(qij.begin(), qij.end());
    return std::move(qst);
  }
  //! Adding other Qshapes
  /*! This is related to direct sum of QSTArray
   *  e.g. { q1, q2 } + { q3, q4 } = { q1, q2, q3, q4 } */
  Qshapes operator+ (const Qshapes& other) const {
    Qshapes qa(*this);
    qa.insert(qa.end(), other.begin(), other.end());
    return std::move(qa);
  }
  //! Adding other Qshapes to this
  Qshapes& operator+= (const Qshapes& other) {
    this->insert(this->end(), other.begin(), other.end());
    return *this;
  }
  //! Return copy of this
  Qshapes operator+ () const {
    Qshapes qp(*this);
    return std::move(qp);
  }
  //! Return conjugated quantum numbers
  Qshapes operator- () const {
    Qshapes qm;
    qm.reserve(this->size());
    for(const Q& qi : *this) qm.push_back(-qi);
    return std::move(qm);
  }
};

//! Indexed product of TVector<Qshapes<Q>, N>
/*! e.g. index = { i, j, k, l }
 *  return vec[0][i] * vec[1][j] * vec[2][k] * vec[3][l] */
template<size_t N, class Q>
inline Q operator* (const TVector<Qshapes<Q>, N>& vec, const IVector<N>& index) {
  Q prod = vec[0][index[0]];
  for(int i = 1; i < N; ++i) prod = prod * vec[i][index[i]];
  return prod;
}

//! Indexed product of TVector<Qshapes<Q>, N>
/*! e.g. index = { i, j, k, l }
 *  return vec[0][i] * vec[1][j] * vec[2][k] * vec[3][l] */
template<size_t N, class Q>
inline Q operator* (const IVector<N>& index, const TVector<Qshapes<Q>, N>& vec) {
  Q prod = vec[0][index[0]];
  for(int i = 1; i < N; ++i) prod = prod * vec[i][index[i]];
  return prod;
}

};

#endif // _BTAS_CXX11_QSHAPES_H
