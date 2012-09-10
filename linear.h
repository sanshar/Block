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

This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_LINEAR_HEADER_H
#define SPIN_LINEAR_HEADER_H
#define WANT_MATH
#define WANT_STREAM
#include <newmat.h>
#include <newmatap.h>
#include <newmatio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include "wavefunction.h"
#include "davidson.h"

namespace SpinAdapted{
namespace Linear
{
  void precondition(Wavefunction& op, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
  void olsenPrecondition(Wavefunction& op, Wavefunction& C0, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
  void block_davidson(vector<Wavefunction>& b, DiagonalMatrix& e, double normtol, const bool &warmUp, Davidson_functor& h_mult, bool& useprecond, bool& solved);
};
}
#endif

