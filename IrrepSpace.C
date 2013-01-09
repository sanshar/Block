/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "IrrepSpace.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace SpinAdapted{

  std::vector<IrrepSpace> IrrepSpace::operator+=(IrrepSpace rhs)
  {
    std::vector<int> irreps = Symmetry::add(irrep, rhs.irrep);
    std::vector<IrrepSpace> vec;
    for (int i=0; i<irreps.size(); i++) {
      vec.push_back(IrrepSpace(irreps[i]));
    }
    return vec;
  }

  bool IrrepSpace::operator==(IrrepSpace rhs) const
  {
    return irrep == rhs.irrep;
  }
  
  bool IrrepSpace::operator!=(IrrepSpace rhs) const
  {
    return irrep != rhs.irrep;
  }
  
  std::vector<IrrepSpace> operator+(IrrepSpace lhs, IrrepSpace rhs)
  {
    return lhs+=rhs;
  }
  
  IrrepSpace operator-(const IrrepSpace lhs)
  {
    if (sym == "trans") {
      int outirrep = Symmetry::negativeof(lhs.irrep);
      return IrrepSpace(outirrep);
    }
    else if (sym == "lzsym")
      return IrrepSpace(-lhs.irrep);
    else
      return lhs;
  }
  
  void IrrepSpace::Save(std::ofstream &ofs) 
  {
    boost::archive::binary_oarchive save_sym(ofs);
    save_sym << *this;
  }
  
  void IrrepSpace::Load (std::ifstream &ifs)
  {
    boost::archive::binary_iarchive load_sym(ifs);
    load_sym >> *this;
  }
  
  ostream& operator<< (ostream& os, const IrrepSpace s)
  {
    if (sym == "dinfh") {
      char goru = s.irrep%2 == 0 ? 'g' : 'u';
      os<< max(0,(s.irrep-2)/2)<<goru;
      if (s.irrep <2) os<< '+';
      else if (s.irrep >=2 && s.irrep <4 ) os<< '-';
    }
    else if (sym == "trans") {
      std::vector<int> irreps = Symmetry::decompress(s.irrep);
      os<<irreps[2]<<irreps[1]<<irreps[0];
    }
    else if (sym == "lzsym")
      os << s.irrep;
    else {
      os<< s.irrep+1;
    }
    return os;
  }
  
  bool IrrepSpace::operator< (const IrrepSpace s) const
  {
    return irrep<s.irrep;
  }

  int IrrepSpace::getirrep() const {return irrep;}
  
}
