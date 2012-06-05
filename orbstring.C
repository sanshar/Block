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


#include "orbstring.h"
#include "global.h"


namespace SpinAdapted{
  namespace OrbstringParams
  {
    int orb_size;
    bool initialised;
  }

}

int SpinAdapted::Orbstring::size () const
{
  return OrbstringParams::orb_size;
}

void SpinAdapted::Orbstring::init (const int o)
{
  OrbstringParams::orb_size = o;
  OrbstringParams::initialised = true;
}

SpinAdapted::Orbstring::Orbstring () {
  sign =0; 
  empty = true;
  occ_rep = std::vector<bool>(OrbstringParams::orb_size, 0); 
}


SpinAdapted::Orbstring& SpinAdapted::Orbstring::c (int i)
  // applies creation operator i to ket
{
  assert ((i >= 0) && (i < OrbstringParams::orb_size));
  if (empty)
    return *this;
  else if (occ_rep[i])
    {
      empty = true;
      return *this;
    }
  else 
    {
      sign *= this->parity(i);
      occ_rep[i] = 1;
    }
  return *this;
}

SpinAdapted::Orbstring& SpinAdapted::Orbstring::d (int i)  // applies destruction operator i to ket
{
  assert ((i >= 0) && (i < OrbstringParams::orb_size));
  if (empty)
    return *this;
  else if (!occ_rep[i])
    {
      empty = true;
      return *this;
    }
  else
    {
      sign *= this->parity(i);
      occ_rep[i] = 0;
    }
  return *this;
}

