/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
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
  void InitStartingBlock (SpinBlock& startingBlock, const bool &forward, int leftState, int rightState,
                          const int & forward_starting_size, const int &backward_starting_size,
                          const int& restartSize, const bool &restart, const bool& warmUp, int integralIndex, const vector<SpinQuantum>& braquanta = vector<SpinQuantum>() , const vector<SpinQuantum>& ketquanta = vector<SpinQuantum>() );

  void InitNewSystemBlock(SpinBlock &system, SpinBlock &systemDot, SpinBlock &newSystem, int leftState, 
			  int rightState, const int &sys_add, const bool &direct, int integralIndex, 
			  const Storagetype &storage, bool haveNormops, bool haveCompops, int constraint=NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

void InitNewSystemBlock(SpinBlock &system, SpinBlock &systemDot, SpinBlock &newSystem, int leftState, int rightState, const int& sys_add, const bool &direct, int integralIndex, const Storagetype &storage, bool haveNormops, bool haveCompops, int constraint, const std::vector<SpinQuantum>& braquanta, const std::vector<SpinQuantum>& ketquanta);

  void InitNewEnvironmentBlock(SpinBlock &environment, SpinBlock& environmentDot, SpinBlock &newEnvironment,
                               const SpinBlock &system, SpinBlock &systemDot, int leftState, int rightState,
                               const int &sys_add, const int &env_add, const bool &forward, const bool &direct, const bool &onedot, 
			       const bool &nexact, const bool &useSlater, int integralIndex, 
			       bool haveNormops, bool haveCompops, const bool& dot_with_sys, int constraint=NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

void InitNewEnvironmentBlock(SpinBlock &environment, SpinBlock& environmentDot, SpinBlock &newEnvironment, 
						      const SpinBlock &system, SpinBlock &systemDot, int leftState, int rightState,
						      const int &sys_add, const int &env_add, const bool &forward, const bool &direct, 
						      const bool &onedot, const bool &nexact, const bool &useSlater, int integralIndex, 
						      bool haveNormops, bool haveCompops, const bool& dot_with_sys, int constraint, const std::vector<SpinQuantum>& braquanta, const std::vector<SpinQuantum>& ketquanta);

  void InitNewOverlapEnvironmentBlock(SpinBlock &environment, SpinBlock& environmentDot, SpinBlock &newEnvironment, 
				      const SpinBlock &system, SpinBlock &systemDot, int leftState, int rightState,
				      const int &sys_add, const int &env_add, const bool &forward, int integralIndex,
				      const bool &onedot, const bool& dot_with_sys, int constraint=NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

  void InitBigBlock(SpinBlock &leftBlock, SpinBlock &rightBlock, SpinBlock &big, const vector<SpinQuantum>& braquanta = vector<SpinQuantum>() , const vector<SpinQuantum>& ketquanta = vector<SpinQuantum>() );

}
}
#endif

