/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SPACE_HEADER
#define SPIN_SPACE_HEADER
#include <boost/serialization/serialization.hpp>
#include <vector>
#include <iostream>
//In nonabelian symmetries, in this case the spin SU(2), each irrep could have many basis vectors, 
//this class represents the irrep (in this case Spin). Each basis vector in a irrep is represented by a 
// different class IrrepVector (Sz). If we dont want to use SpinAdapted algorithm then this class
//will become Sz and each irrep will only have a single basis vector.
namespace SpinAdapted{
class SpinSpace
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & irrep;
  }
  int irrep;
 public:
  SpinSpace() : irrep(0) {}
  explicit SpinSpace(int ir);

  std::vector<SpinSpace> operator+=(SpinSpace rhs);
  bool operator==(SpinSpace rhs) const;
  bool operator!=(SpinSpace rhs) const;
  bool operator<(SpinSpace rhs) const;
  friend std::vector<SpinSpace> operator+(SpinSpace lhs, SpinSpace rhs);
  friend std::vector<SpinSpace> operator-(SpinSpace lhs, SpinSpace rhs);
  friend SpinSpace operator-(SpinSpace lhs);
  void Save(std::ofstream &ofs);
  void Load(std::ifstream &ifs);
  friend std::ostream& operator<<(std::ostream& os, const SpinSpace s);
  int getirrep() const ;
};
}
#endif
