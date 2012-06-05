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

#ifndef SPIN_DAVIDSON_HEADER
#define SPIN_DAVIDSON_HEADER
#include "wavefunction.h"
#include "spinblock.h"

namespace SpinAdapted{
struct Davidson_functor
{
  virtual void operator()(Wavefunction& c, Wavefunction& v) = 0;
  virtual const SpinBlock& get_block() = 0;
};

class multiply_h : public Davidson_functor
{
 private:
  const SpinBlock& block;
 public:
  multiply_h(const SpinBlock& b, const bool &onedot_);
  void operator()(Wavefunction& c, Wavefunction& v);
  const SpinBlock& get_block() {return block;}
};
}

#endif
