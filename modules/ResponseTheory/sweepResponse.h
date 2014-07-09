/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SWEEPRESPONSE_HEADER
#define SPIN_SWEEPRESPONSE_HEADER
#include "spinblock.h"
#include "sweep_params.h"


namespace SpinAdapted{
namespace SweepResponse
{
  void BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, int correctionVector, int baseState);
  double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int targetState, int correctionVector, int baseState);
  void StartUp (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool& dot_with_sys, int targetState, int correctionVector, int baseState);
  void WavefunctionCanonicalize (SweepParams &sweepParams, SpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int targetState, int correctionVector, int baseState);
}
}
#endif

