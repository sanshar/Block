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


#ifndef SPIN_SOLVER_HEADER_H
#define SPIN_SOLVER_HEADER_H
#include "spinblock.h"
#include "wavefunction.h"

namespace SpinAdapted{
namespace Solver
{
  void solve_wavefunction(vector<Wavefunction>& solution, std::vector<double>& energies, SpinBlock& big, const double tol, 
			  const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp, double additional_noise=0.0);
};
}
#endif
