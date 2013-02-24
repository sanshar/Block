/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SWEEP_LRT_HEADER
#define SPIN_SWEEP_LRT_HEADER
#include "module/dmrg_lrt/spinblock_deriv.h"
#include "sweep_params.h"

namespace SpinAdapted{

namespace LRT
{

void BlockAndDecimate(SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys);

double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);

};

};
#endif

