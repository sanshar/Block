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
#ifndef DEC
using namespace std;
#endif

namespace SpinAdapted{
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
  int totalSpin;
  IrrepSpace orbitalSymmetry;

  SpinQuantum ();
  SpinQuantum (const int p, const int s, const SpinAdapted::IrrepSpace orbS);
  SpinQuantum operator-() const;
  //SpinQuantum operator- (const SpinQuantum q) const;
  int insertionNum(const SpinQuantum& ql, const SpinQuantum& qr) const;
  vector<SpinQuantum> operator+ (const SpinQuantum q) const;
  vector<SpinQuantum> operator- (const SpinQuantum q) const;
  bool operator== (const SpinQuantum q) const;
  bool operator!= (const SpinQuantum q) const;
  bool operator< (const SpinQuantum q) const;
  friend bool IsFermion (const SpinQuantum q);
  friend ostream& operator<< (ostream& os, const SpinQuantum q);
  void Save (std::ofstream &ofs);
  void Load (std::ifstream &ifs);
  int get_s() const {return totalSpin;}
  int get_n() const {return particleNumber;}
  SpinAdapted::IrrepSpace get_symm() const {return orbitalSymmetry;}
  bool allow(const SpinQuantum s1, const SpinQuantum s2) const;

  static bool can_complement (SpinQuantum q);
  void complementize ();
  vector<SpinQuantum> get_complement() const;
};
}  


#endif

