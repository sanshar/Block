
//
/*! \file  Quantum.h
 *  \brief Default quantum number class named Quantum.
 */

#ifndef _BTAS_CXX11_DEFAULT_QUANTUM_H
#define _BTAS_CXX11_DEFAULT_QUANTUM_H 1

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>

//! Default quantum number class
/*! Minimum design of quantum number class provided in BTAS (designed as a spin quantum number).
 *
 *  Quantum number class requires,
 *  - Boost serialization function
 *  - Global function which returns zero quantum number
 *  - Default constructor
 *  - Comparison operators, ==, !=, <, and >
 *  - Multiplication operator * to take product of two quantum number
 *  - Function which returns parity (true for scaling by -1)
 *  - Function which returns Clebsch-Gordan coefficient
 *  - Unary plus and minus operators for conjugation
 *  - Stream printing operator<<
 *
 *  If defined _DEFAULT_QUANTUM, this is included as a default quantum number class.
 */
class Quantum {
private:
  friend class boost::serialization::access;
  //! Boost serialization
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_q_val; }
public:
  //! Global function to return zero quantum number
  const static Quantum zero() { return Quantum(0); }
  //! Default constructor
  Quantum() : m_q_val(0) { }
  //! Initializer
  Quantum(int q_val) : m_q_val(q_val) { }
  //! Equal to 
  inline bool operator== (const Quantum& other) const { return m_q_val == other.m_q_val; }
  //! Not equal to
  inline bool operator!= (const Quantum& other) const { return m_q_val != other.m_q_val; }
  //! Less than
  inline bool operator<  (const Quantum& other) const { return m_q_val <  other.m_q_val; }
  //! Greater than
  inline bool operator>  (const Quantum& other) const { return m_q_val >  other.m_q_val; }
  //! Product of two quantum number
  inline Quantum operator* (const Quantum& other) const { return Quantum(m_q_val+other.m_q_val); }
  //! Parity, i.e. return -1 if fermion and odd particle number
  inline bool   parity()  const { return false; }
  //! Clebsch-Gordan coefficient for non-Abelian symmetry
  inline double clebsch() const { return 1.0; }
  //! Unary plus
  Quantum operator+ () const { return Quantum(+m_q_val); }
  //! Unary minus to return conjugated quantum number
  Quantum operator- () const { return Quantum(-m_q_val); }
  //! Printing function
  friend std::ostream& operator<< (std::ostream& ost, const Quantum& q) {
    ost << "(" << std::setw(2) << q.m_q_val << ")";
    return ost;
  }
private:
  //! Quantum number
  int m_q_val;
};

#endif // _BTAS_CXX11_DEFAULT_QUANTUM_H
