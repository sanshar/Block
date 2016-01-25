/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPINQUANTUM_HEADER
#define SPINQUANTUM_HEADER
#include <iostream>
#include <stdio.h>
#include "Symmetry.h"
#include "IrrepSpace.h"
#include <vector>
#include <boost/serialization/serialization.hpp>
#include "SpinSpace.h"
#ifndef DEC
using namespace std;
#endif


namespace SpinAdapted{

/// @brief Quantum numbers in DMRG
///
/// SpinQuantum carries particle number, orbital symmetry (Abelian and non-Abelian)
/// and Spin (Sz or S2) labels. Arithmetic (addition, subtraction, comparison) of
/// quantum numbers is defined. For non-Abelian groups, addition of two quantum
/// numbers produces a set of quantum numbers (e.g. a Clebsch-Gordon series).
class SpinQuantum 
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & particleNumber & totalSpin & orbitalSymmetry;
  }
  
 public:
  int particleNumber;
  SpinSpace totalSpin;
  IrrepSpace orbitalSymmetry;

  SpinQuantum ();
  SpinQuantum (const int p, const SpinSpace s, const SpinAdapted::IrrepSpace orbS);
  SpinQuantum operator-() const;

  /// \return Vector of (Sz type) quantum numbers with Sz=-S ... +S
  vector<SpinQuantum> spinToNonSpin() const;
  int insertionNum(const SpinQuantum& ql, const SpinQuantum& qr) const;

  /// Adds quantum number (to self)
  /// \return Clebsch-Gordon sum of quantum numbers (if non-Abelian)
  /// or length-1 vector if Abelian.
  vector<SpinQuantum> operator+ (const SpinQuantum& q) const;
  vector<SpinQuantum> operator+ (const vector<SpinQuantum>& q) const;

  vector<SpinQuantum> operator- (const SpinQuantum& q) const;
  vector<SpinQuantum> operator- (const vector<SpinQuantum>& q) const;

  bool operator== (const SpinQuantum q) const;
  bool operator!= (const SpinQuantum q) const;
  bool operator< (const SpinQuantum q) const;

  /// \return True if particleNumber is odd..
  friend bool IsFermion (const SpinQuantum q);
  friend ostream& operator<< (ostream& os, const SpinQuantum q);
  void Save (std::ofstream &ofs);
  void Load (std::ifstream &ifs);
  SpinAdapted::SpinSpace get_s() const {return totalSpin;}
  int get_n() const {return particleNumber;}
  SpinAdapted::IrrepSpace get_symm() const {return orbitalSymmetry;}

  /// \return true if self is contained in (Clebsch-Gordon) sum, s1+s2.
  bool allow(const SpinQuantum s1, const SpinQuantum s2) const;


  static bool can_complement (SpinQuantum q);
  vector<SpinQuantum> get_complement() const;
};
}  


#endif

