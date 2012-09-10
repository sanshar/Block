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


#ifndef SPIN_INIT_BLOCKS_HEADER
#define SPIN_INIT_BLOCKS_HEADER
#include "spinblock.h"
#include "sweep_params.h"

namespace SpinAdapted{

namespace InitBlocks
{
  void InitStartingBlock (SpinBlock& startingBlock, const bool &forward,
                          const int & forward_starting_size, const int &backward_starting_size,
                          const int& restartSize, const bool &restart, const bool& warmUp);
  void InitNewSystemBlock(SpinBlock &system, SpinBlock &systemDot, SpinBlock &newSystem, const int &sys_add, const bool &direct, const Storagetype &storage= DISTRIBUTED_STORAGE, bool haveNormops = true, bool haveCompops = true);

  void InitNewEnvironmentBlock(SpinBlock &environment, SpinBlock& environmentDot, SpinBlock &newEnvironment,
                               const SpinBlock &system, SpinBlock &systemDot,
                               const int &sys_add, const int &env_add, const bool &forward, const bool &direct, const bool &onedot, 
			       const bool &nexact, const bool &useSlater, 
			       bool haveNormops = true, bool haveCompops = true, const bool& dot_with_sys = true);
  void InitBigBlock(SpinBlock &leftBlock, SpinBlock &rightBlock, SpinBlock &big);
}
}
#endif

