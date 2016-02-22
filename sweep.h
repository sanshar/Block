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
  void do_overlap(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);
  void fullci(double sweep_tol);
  void tiny(double sweep_tol);

  void CanonicalizeWavefunction(SweepParams &sweepParams, const bool &forward, int currentstate);
  void InitializeStateInfo(SweepParams &sweepParams, const bool &forward, int currentstate);
  void InitializeOverlapSpinBlocks(SweepParams &sweepParams, const bool &forward, int stateA, int stateB, int integralIndex);
  void calculateAllOverlap(Matrix& overlap);
  void calculateHMatrixElements(Matrix& H);
  void makeSystemEnvironmentBigBlocks(SpinBlock& system, SpinBlock& systemDot, SpinBlock& newSystem, 
				      SpinBlock& environment, SpinBlock& environmentDot, SpinBlock& newEnvironment,
				      SpinBlock& big, SweepParams& sweepParams, const bool& dot_with_sys, const bool& useSlater,
				      int integralIndex, int braState=-1, int ketState=-1, const vector<SpinQuantum>& braquanta= vector<SpinQuantum>(), const vector<SpinQuantum>& ketquanta= vector<SpinQuantum>());
  void makeSystemEnvironmentBigOverlapBlocks(const std::vector<int>& systemSites, SpinBlock& systemDot, SpinBlock& environmentDot,
					     SpinBlock& system, SpinBlock& newSystem, SpinBlock& environment, SpinBlock& newEnvironment,
					     SpinBlock& big, SweepParams& sweepParams, const bool& dot_with_sys, const bool& useSlater,
					     int integralIndex, int braState, int ketState);
  
};
}
#endif

