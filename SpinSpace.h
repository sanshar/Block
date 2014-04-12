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


namespace SpinAdapted{
/// @brief `SpinSpace` is a wrapper class for the spin irrep.
///
/// In Sz symmetry, the irrep is 2*Sz. In SU(2) symmetry, the
/// irrep is 2*S. The behaviour is toggled between Sz and S2 symmetry by checking dmrg.spinAdapted().
///
/// `SpinSpace` defines addition of irreps in SU(2) symmetry
/// to return the Clebsch-Gordon series for S1+S2=|S1-S2|...S1+S2.
///
/// 
class SpinSpace
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & irrep;
  }
  /// Integer representing 2*S or 2*Sz 
  int irrep;
 public:
  SpinSpace() : irrep(0) {}
  explicit SpinSpace(int ir);

  std::vector<SpinSpace> operator+=(SpinSpace rhs);
  bool operator==(SpinSpace rhs) const;
  bool operator!=(SpinSpace rhs) const;
  bool operator<(SpinSpace rhs) const;
  /// Adds integer irreps in lhs, rhs. 
  /// \return If S2 symmetry (`dmrg.spinAdapted()==true`),  vector |S1-S2| ... S1+S2, else, 
  /// vector of length 1 containing Sz1+Sz2
  friend std::vector<SpinSpace> operator+(SpinSpace lhs, SpinSpace rhs);

  /// Negate irrep.
  /// \return In Sz symmetry, -Sz; for S2 symmetry, returns same irrep S (i.e. does nothing,
  /// since negative S has no meaning).
  friend SpinSpace operator-(SpinSpace lhs);
  void Save(std::ofstream &ofs);
  void Load(std::ifstream &ifs);
  friend std::ostream& operator<<(std::ostream& os, const SpinSpace s);
  int getirrep() const ;
};
}
#endif
