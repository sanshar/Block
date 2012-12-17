/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "global.h"
#include "initblocks.h"
#include "pario.h"

void SpinAdapted::InitBlocks::InitStartingBlock (SpinBlock& startingBlock, const bool &forward, 
                                    const int & forward_starting_size, const int &backward_starting_size,
                                    const int& restartSize, const bool &restart, const bool& warmUp)
{
  if (restart && restartSize != 1)
  {
    int len = restart? restartSize : forward_starting_size;
    vector<int> sites(len);
    if (forward)
      for (int i=0; i<len; i++)
	sites[i] = i;
    else
      for (int i=0; i<len; i++) 
	sites[i] = dmrginp.last_site() - len +i ;
    
    if (restart)
      SpinBlock::restore (forward, sites, startingBlock);
    else
      SpinBlock::restore (true, sites, startingBlock);
  }
  else if (forward)
  {
    startingBlock = SpinBlock(0, forward_starting_size - 1, true);
    
    if (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s() != 0)
    {
      SpinQuantum s = dmrginp.molecule_quantum();
      s = SpinQuantum(s.get_s(), s.get_s(), IrrepSpace(0));
      int qs = 1, ns = 1;
      StateInfo addstate(ns, &s, &qs); 
      SpinBlock dummyblock(addstate);
      SpinBlock newstartingBlock;
      newstartingBlock.default_op_components(false, startingBlock, dummyblock, true, true);
      newstartingBlock.setstoragetype(LOCAL_STORAGE);
      newstartingBlock.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, startingBlock, dummyblock);
      startingBlock.clear();
      startingBlock = newstartingBlock;
    }
  }
  else
  {
    std::vector<int> backwardSites;
    for (int i = 0; i < backward_starting_size; ++i) 
	    backwardSites.push_back (dmrginp.last_site() - i - 1);
    sort (backwardSites.begin (), backwardSites.end ());
	  startingBlock.default_op_components(false);
    startingBlock.BuildTensorProductBlock (backwardSites);
  }
}


void SpinAdapted::InitBlocks::InitNewSystemBlock(SpinBlock &system, SpinBlock &systemDot, SpinBlock &newSystem, const int &sys_add, const bool &direct, const Storagetype &storage, bool haveNormops, bool haveCompops)
{

  newSystem.default_op_components(direct, system, systemDot, haveNormops, haveCompops);
  newSystem.setstoragetype(storage);
  newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemDot);

  if (dmrginp.outputlevel() > 0) {
    pout << "\t\t\t NewSystem block " << endl << newSystem << endl;
    newSystem.printOperatorSummary();
  }
}

void SpinAdapted::InitBlocks::InitNewEnvironmentBlock(SpinBlock &environment, SpinBlock& environmentDot, SpinBlock &newEnvironment, 
                                          const SpinBlock &system, SpinBlock &systemDot,
					 const int &sys_add, const int &env_add, const bool &forward, const bool &direct, const bool &onedot, const bool &nexact, const bool &useSlater, bool haveNormops, bool haveCompops, const bool& dot_with_sys)
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
    StateInfo system_stateinfo = system.get_stateInfo();
    StateInfo sysdot_stateinfo = systemDot.get_stateInfo();
    StateInfo tmp;
    TensorProduct (system_stateinfo, sysdot_stateinfo, tmp, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    tmp.CollectQuanta ();

    if (dmrginp.do_fci() || environmentSites.size() == nexact)
    {
      if ((!dot_with_sys && onedot) || !onedot)
      {
	environment.default_op_components(!forward);
	environment.setstoragetype(DISTRIBUTED_STORAGE);
	environment.BuildTensorProductBlock(environmentSites);
	SpinBlock::store (true, environmentSites, environment);	
      }
      else
      {
	newEnvironment.default_op_components(!forward);
	newEnvironment.setstoragetype(DISTRIBUTED_STORAGE);
	newEnvironment.BuildTensorProductBlock(environmentSites);
	SpinBlock::store (true, environmentSites, newEnvironment);	
      }

    }
    else //used for warmup guess environemnt
    {
      std::vector<SpinQuantum> quantumNumbers;
      std::vector<int> distribution;
      std::map<SpinQuantum, int> quantaDist;
      std::map<SpinQuantum, int>::iterator quantaIterator;
      bool environmentComplementary = !forward;
      StateInfo tmp2;


      if (onedot)// || useSlater)
	tmp.quanta_distribution (quantumNumbers, distribution, true);
      else 
      {
        StateInfo environmentdot_stateinfo = environmentDot.get_stateInfo();
        TensorProduct (tmp, environmentdot_stateinfo, tmp2, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
        tmp2.CollectQuanta ();
        tmp2.quanta_distribution (quantumNumbers, distribution, true);
	//tmp.quanta_distribution (quantumNumbers, distribution, true);
      }

      for (int i = 0; i < distribution.size (); ++i)
      {
	quantaIterator = quantaDist.find(quantumNumbers[i]);
	if (quantaIterator != quantaDist.end())
	  distribution[i] += quantaIterator->second;

        distribution [i] /= 4; distribution [i] += 1;
        if (distribution [i] > dmrginp.nquanta()){
	  distribution [i] = dmrginp.nquanta();
	}
	
	if(quantaIterator != quantaDist.end())
	  quantaIterator->second = distribution[i];
	else
	  quantaDist[quantumNumbers[i]] = distribution[i];
      }
      if (dmrginp.outputlevel() > 0)
	pout << "\t\t\t Quantum numbers and states used for warm up :: " << endl << "\t\t\t ";
      quantumNumbers.clear(); quantumNumbers.reserve(distribution.size());
      distribution.clear();distribution.reserve(quantumNumbers.size());
      std::map<SpinQuantum, int>::iterator qit = quantaDist.begin();

      for (; qit != quantaDist.end(); qit++)
      {
	quantumNumbers.push_back( qit->first); distribution.push_back(qit->second); 
	if (dmrginp.outputlevel() > 0) {
	  pout << quantumNumbers.back() << " = " << distribution.back() << ", ";
	  if (! (quantumNumbers.size() - 6) % 6) pout << endl << "\t\t\t ";
	}
      }
      pout << endl;


      if(dot_with_sys && onedot)
      {
        newEnvironment.BuildSlaterBlock (environmentSites, quantumNumbers, distribution, false, false);
      }
      else 
      {
        environment.BuildSlaterBlock (environmentSites, quantumNumbers, distribution, false, haveNormops);
      }
    }
  }
  else
  {
    if (dmrginp.outputlevel() > 0)
      pout << "\t\t\t Restoring block of size " << environmentSites.size () << " from previous iteration" << endl;
    if(dot_with_sys && onedot)
      SpinBlock::restore (!forward, environmentSites, newEnvironment);
    else
      SpinBlock::restore (!forward, environmentSites, environment);
    if (dmrginp.outputlevel() > 0)
      mcheck("");
  }

  // now initialise newEnvironment
  if (!dot_with_sys || !onedot)
  {
    dmrginp.datatransfer -> start();
    environment.addAdditionalCompOps();
    dmrginp.datatransfer -> stop();

      newEnvironment.default_op_components(direct, environment, environmentDot, haveNormops, haveCompops);
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

void SpinAdapted::InitBlocks::InitBigBlock(SpinBlock &leftBlock, SpinBlock &rightBlock, SpinBlock &big)
{
  //set big block components
  big.set_big_components(); 
  // build the big block
  big.BuildSumBlock(PARTICLE_SPIN_NUMBER_CONSTRAINT, leftBlock, rightBlock);
}
