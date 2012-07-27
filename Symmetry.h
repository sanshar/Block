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


#ifndef SPIN_SYMMETRY_HEADER
#define SPIN_SYMMETRY_HEADER
#include <fstream>
#include <cstdio>
#include <multiarray.h>
#include <boost/serialization/serialization.hpp>
#include <vector>
#include "IrrepVector.h"
#include <string.h>

using namespace boost;
/** Non-Abelian and Abelian Symmetry class */

namespace SpinAdapted {
class IrrepVector;

class Symmetry
{
 public:
  static void InitialiseTable(string sym);

  static std::vector<int> add(int irrepl, int irrepr);
  static int sizeofIrrep(int irrep);
  static bool irrepAllowed(int irrep);
  static string stringOfIrrep(int irrep);
  static double spatial_cg(int a, int b, int c, int la, int lb, int lc) ;
  static double spatial_sixj(int a, int b, int c, int d, int e, int f) ;
  static double spatial_ninej(int j1, int j2, int j12, int j3, int j4, int j34, int j13, int j24, int j);
  static std::vector<int> decompress(int irrep);
  static int compress(std::vector<int>& irreps);
};
}
#endif
