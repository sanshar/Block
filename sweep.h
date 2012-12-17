/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SWEEP_HEADER
#define SPIN_SWEEP_HEADER
#include "spinblock.h"
#include "sweep_params.h"

namespace SpinAdapted{
namespace Sweep
{
  void BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys);
  void Startup (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem);
  double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);
  void fullci(double sweep_tol);
  void tiny(double sweep_tol);
};
}
#endif

