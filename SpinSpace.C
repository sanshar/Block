/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "SpinSpace.h"
#include "global.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace SpinAdapted{

  SpinSpace::SpinSpace(int ir) {
    if (dmrginp.spinAdapted())
      irrep = abs(ir);
    else
      irrep = ir;
  }

  std::vector<SpinSpace> SpinSpace::operator+=(SpinSpace rhs)
  {
    std::vector<SpinSpace> Spins;
    if (!dmrginp.spinAdapted()) 
      Spins.push_back(SpinSpace(irrep+rhs.irrep));
    else {
      for (int i=abs(irrep-rhs.irrep); i<=irrep+rhs.irrep; i+=2)
	Spins.push_back(SpinSpace(i));
    }

    return Spins;
  }

  bool SpinSpace::operator==(SpinSpace rhs) const
  {
    return irrep == rhs.irrep;
  }
  
  bool SpinSpace::operator!=(SpinSpace rhs) const
  {
    return irrep != rhs.irrep;
  }
  
  std::vector<SpinSpace> operator+(SpinSpace lhs, SpinSpace rhs)
  {
    return lhs+=rhs;
  }
  
  SpinSpace operator-(const SpinSpace lhs)
  {
    int outirrep;
    if (!dmrginp.spinAdapted()) 
      outirrep = -lhs.irrep;
    else
      outirrep = lhs.irrep;
    return SpinSpace(outirrep);
  }
  
  void SpinSpace::Save(std::ofstream &ofs) 
  {
    boost::archive::binary_oarchive save_sym(ofs);
    save_sym << *this;
  }
  
  void SpinSpace::Load (std::ifstream &ifs)
  {
    boost::archive::binary_iarchive load_sym(ifs);
    load_sym >> *this;
  }
  
  ostream& operator<< (ostream& os, const SpinSpace s)
  {
    os << s.irrep;
    return os;
  }
  
  bool SpinSpace::operator< (const SpinSpace s) const
  {
    return irrep<s.irrep;
  }

  int SpinSpace::getirrep() const {return irrep;}
  
}
