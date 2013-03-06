/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef LRT_SPIN_INIT_BLOCKS_HEADER
#define LRT_SPIN_INIT_BLOCKS_HEADER
#include "spinblock.h"
#include "sweep_params.h"

namespace SpinAdapted {

namespace InitBlocks {

namespace LRT {

void InitNewEnvironmentBlock
(SpinBlock &environment, SpinBlock& environmentDot, SpinBlock &newEnvironment, const SpinBlock &system, SpinBlock &systemDot,
 const Matrix& alpha, const int &sys_add, const int &env_add, const bool &forward, const bool &direct, const bool &onedot, 
 const bool &nexact, const bool &useSlater, bool haveNormops = true, bool haveCompops = true, const bool& dot_with_sys = true, int nroots = 1);

};

};

};

#endif

