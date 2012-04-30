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
    else {
      os<< s.irrep;
    }
    return os;
  }
  
  bool IrrepSpace::operator< (const IrrepSpace s) const
  {
    return irrep<s.irrep;
  }

  int IrrepSpace::getirrep() const {return irrep;}
  
}
