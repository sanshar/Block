/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "global.h"
#include "modules/lrt/lrt_initblocks.h"
#include "pario.h"

void SpinAdapted::InitBlocks::LRT::InitNewEnvironmentBlock
(SpinBlock& environment, SpinBlock& environmentDot, SpinBlock& newEnvironment, const SpinBlock& system, SpinBlock& systemDot,
 const Matrix& alpha, const int& sys_add, const int& env_add, const bool& forward, const bool& direct, const bool& onedot,
 const bool& nexact, const bool& useSlater, bool haveNormops, bool haveCompops, const bool& dot_with_sys, int nroots)
{
  // now initialise environment Dot
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int systemDotSize = sys_add - 1;
  int environmentDotSize = env_add - 1;
  if (forward)
  {
    systemDotStart = *system.get_sites().rbegin () + 1;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
    environmentStart = environmentDotEnd + 1;
    environmentEnd = dmrginp.last_site() - 1;
  }
  else
  {
    systemDotStart = system.get_sites() [0] - 1;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
    environmentStart = environmentDotEnd - 1;
    environmentEnd = 0;
  }

  std::vector<int> environmentSites;
  environmentSites.resize(abs(environmentEnd - environmentStart) + 1);
  for (int i = 0; i < abs(environmentEnd - environmentStart) + 1; ++i) *(environmentSites.begin () + i) = min(environmentStart,environmentEnd) + i;


  // now initialise environment
  if (useSlater)
  {
    pout << "\t\t\t InitBlocks::LRT::InitNewEnvironmentBlock  invalid option was specified, could not construct block from slater det." << endl;
    abort();
  }
  else
  {
    if (dmrginp.outputlevel() > 0)
      pout << "\t\t\t Restoring block of size " << environmentSites.size () << " from previous iteration" << endl;
    if(dot_with_sys && onedot) {
      SpinBlock::restore (!forward, environmentSites, newEnvironment);
      newEnvironment.rotatebyRitzVectors(alpha, nroots);
    }
    else {
      SpinBlock::restore (!forward, environmentSites, environment);
      environment.rotatebyRitzVectors(alpha, nroots);
    }
    if (dmrginp.outputlevel() > 0)
      mcheck("");
  }

  // now initialise newEnvironment
  if (!dot_with_sys || !onedot)
  {
    dmrginp.datatransfer -> start();
    // broadcast additional compops (0-th)
    environment.addAdditionalCompOps();
    // broadcast additional compops (1-st)
    for(int i = 1; i < nroots; ++i) {
      environment.addAdditionalCompOps(0, i);
      environment.addAdditionalCompOps(i, 0);
    }
    dmrginp.datatransfer -> stop();

    newEnvironment.default_op_components(direct, environment, environmentDot, haveNormops, haveCompops, nroots);
    newEnvironment.setstoragetype(DISTRIBUTED_STORAGE);
    
    newEnvironment.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, environment, environmentDot);
    if (dmrginp.outputlevel() > 0) {
      pout << "\t\t\t Environment block " << endl << environment << endl;
      environment.printOperatorSummary();
      pout << "\t\t\t NewEnvironment block " << endl << newEnvironment << endl;
      newEnvironment.printOperatorSummary();
    }
  }
  else  if (dmrginp.outputlevel() > 0) {
    pout << "\t\t\t Environment block " << endl << newEnvironment << endl;
    newEnvironment.printOperatorSummary();
  }

}

