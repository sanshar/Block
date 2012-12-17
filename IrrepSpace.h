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


#ifndef IRREP_SPACE_HEADER
#define IRREP_SPACE_HEADER
#include "Symmetry.h"

//In nonabelian symmetries each irrep could have many basis vectors, this class
// represents the irrep. Each basis vector in a irrep is represented by a 
// different class IrrepVector
namespace SpinAdapted{
class IrrepSpace
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & irrep;
  }
  int irrep;
 public:
  IrrepSpace() : irrep(0) {}
  explicit IrrepSpace(int ir) : irrep(ir) {}

  std::vector<IrrepSpace> operator+=(IrrepSpace rhs);
  bool operator==(IrrepSpace rhs) const;
  bool operator!=(IrrepSpace rhs) const;
  bool operator<(IrrepSpace rhs) const;
  friend std::vector<IrrepSpace> operator+(IrrepSpace lhs, IrrepSpace rhs);
  friend std::vector<IrrepSpace> operator-(IrrepSpace lhs, IrrepSpace rhs);
  friend IrrepSpace operator-(IrrepSpace lhs);
  void Save(std::ofstream &ofs);
  void Load(std::ifstream &ifs);
  friend ostream& operator<<(ostream& os, const IrrepSpace s);
  int getirrep() const ;
};
}
#endif
