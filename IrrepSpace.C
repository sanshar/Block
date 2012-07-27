/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        
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
