/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "IntegralMatrix.h"
#include <fstream>
#include "input.h"
#include "pario.h"
#include "global.h"
#include "orbstring.h"
#include "least_squares.h"
#include <include/communicate.h>
#ifdef _OPENMP
#include <omp.h>
#endif
//the following can be removed later
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/format.hpp>
#include "spinblock.h"
#include "StateInfo.h"
#include "operatorfunctions.h"
#include "wavefunction.h"
#include "solver.h"
#include "davidson.h"
#include "guess_wavefunction.h"
#include "rotationmat.h"
#include "density.h"
#include "sweep.h"
#include "sweepgenblock.h"
#include "sweeponepdm.h"
#include "sweeptwopdm.h"
#include "BaseOperator.h"
#include "dmrg_wrapper.h"

#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "pario.h"


void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
void ReadInput(char* conf);
void fullrestartGenblock();
void license() {
#ifndef MOLPRO
  pout << "Block  Copyright (C) 2012  Garnet K.-L. Chan"<<endl;
  pout << "This program comes with ABSOLUTELY NO WARRANTY; for details see license file."<<endl;
  pout << "This is free software, and you are welcome to redistribute it"<<endl;
  pout << "under certain conditions; see license file for details."<<endl;
#endif
}


namespace SpinAdapted{
  Timer globaltimer(false);
  bool DEBUGWAIT = false;
  bool DEBUG_MEMORY = false;
  bool restartwarm = false;
  OneElectronArray v_1;
  TwoElectronArray v_2(TwoElectronArray::restrictedNonPermSymm);
  Input dmrginp;
  int MAX_THRD = 1;
  bool FULLRESTART;
  bool RESTART;
  bool reset_iter;
  std::vector<int> NPROP;
  int PROPBITLEN=1;
}

using namespace SpinAdapted;

int calldmrg(char* input, char* output)
{
  license();
  if (output != 0) {
    ofstream file;
    file.open(output);
    cout.rdbuf(file.rdbuf());
  }

  ReadInput(input);
  MAX_THRD = dmrginp.thrds_per_node()[mpigetrank()];
#ifdef _OPENMP
  omp_set_num_threads(MAX_THRD);
#endif

   //Initializing timer calls
  dmrginp.initCumulTimer();

  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  SweepParams sweepParams;

  SweepParams sweep_copy;
  bool direction_copy; int restartsize_copy;


  switch(dmrginp.calc_type()) {
    
  case (DMRG):
    if (RESTART && !FULLRESTART)
      restart(sweep_tol, reset_iter);
    else if (FULLRESTART) {
      fullrestartGenblock();
      reset_iter = true;
      sweepParams.restorestate(direction, restartsize);
      sweepParams.calc_niter();
      sweepParams.savestate(direction, restartsize);
      restart(sweep_tol, reset_iter);
    }
    else {
      dmrg(sweep_tol);
    }
    break;

  case (FCI):
    Sweep::fullci(sweep_tol);
    break;
    
  case (TINYCALC):
    Sweep::tiny(sweep_tol);
    break;

  case (ONEPDM):
    if (dmrginp.algorithm_method() == TWODOT) {
      pout << "Onepdm not allowed with twodot algorithm" << endl;
      abort();
    }
    if (RESTART && !FULLRESTART)
      restart(sweep_tol, reset_iter);
    else if (FULLRESTART) {
      fullrestartGenblock();
      reset_iter = true;
      sweepParams.restorestate(direction, restartsize);
      sweepParams.calc_niter();
      sweepParams.savestate(direction, restartsize);
      restart(sweep_tol, reset_iter);
    }
    else {
      dmrg(sweep_tol);
    }

    dmrginp.screen_tol() = 0.0; //need to turn screening off for onepdm
    dmrginp.Sz() = dmrginp.total_spin_number();
    dmrginp.do_cd() = true;
    dmrginp.screen_tol() = 0.0;

    sweep_copy.restorestate(direction_copy, restartsize_copy);
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
    dmrginp.set_fullrestart() = false;    

    SweepOnepdm::do_one(sweepParams, false, direction, false, 0);
    sweep_copy.savestate(direction_copy, restartsize_copy);
    break;

  case (TWOPDM):
    if (dmrginp.algorithm_method() == TWODOT) {
      pout << "Twopdm not allowed with twodot algorithm" << endl;
      abort();
    }

    if (RESTART && !FULLRESTART)
      restart(sweep_tol, reset_iter);
    else if (FULLRESTART) {
      fullrestartGenblock();
      reset_iter = true;
      sweepParams.restorestate(direction, restartsize);
      sweepParams.calc_niter();
      sweepParams.savestate(direction, restartsize);
      restart(sweep_tol, reset_iter);
    }
    else {
      dmrg(sweep_tol);
    }


    dmrginp.screen_tol() = 0.0; //need to turn screening off for onepdm
    dmrginp.Sz() = dmrginp.total_spin_number();
    dmrginp.do_cd() = true;
    dmrginp.screen_tol() = 0.0;
    sweep_copy.restorestate(direction_copy, restartsize_copy);

    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
    dmrginp.set_fullrestart() = false;    

    for (int state=0; state<dmrginp.nroots(); state++) {

      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state);
    }
    sweep_copy.savestate(direction_copy, restartsize_copy);
    break;

  case (RESTART_ONEPDM):
    if(sym == "dinfh") {
      pout << "One pdm not implemented with dinfh symmetry"<<endl;
      abort();
    }
    sweepParams.restorestate(direction, restartsize);
    if(!sweepParams.get_onedot() || dmrginp.algorithm_method() == TWODOT) {
      pout << "Onepdm only runs for the onedot algorithm" << endl;
      abort();
    }


    dmrginp.screen_tol() = 0.0; //need to turn screening off for onepdm
    dmrginp.Sz() = dmrginp.total_spin_number();
    dmrginp.do_cd() = true;
    dmrginp.screen_tol() = 0.0;

    sweep_copy.restorestate(direction_copy, restartsize_copy);
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
    dmrginp.set_fullrestart() = false;    

    SweepOnepdm::do_one(sweepParams, false, direction, false, 0);
    sweep_copy.savestate(direction_copy, restartsize_copy);

    break;

  case (RESTART_TWOPDM):
    if(sym == "dinfh") {
      pout << "Two pdm not implemented with dinfh symmetry"<<endl;
      abort();
    }
    sweepParams.restorestate(direction, restartsize);
    if(!sweepParams.get_onedot() || dmrginp.algorithm_method() == TWODOT) {
      pout << "Twopdm only runs for the onedot algorithm" << endl;
      abort();
    }

    dmrginp.screen_tol() = 0.0; //need to turn screening off for onepdm
    dmrginp.Sz() = dmrginp.total_spin_number();
    dmrginp.do_cd() = true;
    dmrginp.screen_tol() = 0.0;    
    sweep_copy.restorestate(direction_copy, restartsize_copy);

    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
    dmrginp.set_fullrestart() = false;    

    for (int state=0; state<dmrginp.nroots(); state++) {

      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state);
    }
    sweep_copy.savestate(direction_copy, restartsize_copy);

    break;

  // Run into DMRG Linear-Response Theory for Excited State
  case (DMRG_LRT):
    if (dmrginp.algorithm_method() == TWODOT) {
      pout << "DMRG-LRT not allowed with twodot algorithm" << endl;
      abort();
    }

    // TODO: followings are just copied from other calc_type().
    //       figure out all necessary setup before DMRG-LRT
    if (RESTART && !FULLRESTART)
      restart(sweep_tol, reset_iter);
    else if (FULLRESTART) {
      fullrestartGenblock();
      reset_iter = true;
      sweepParams.restorestate(direction, restartsize);
      sweepParams.calc_niter();
      sweepParams.savestate(direction, restartsize);
      restart(sweep_tol, reset_iter);
    }

    dmrginp.screen_tol() = 0.0;
    dmrginp.Sz() = dmrginp.total_spin_number();
    dmrginp.do_cd() = true;
    dmrginp.screen_tol() = 0.0;

    sweep_copy.restorestate(direction_copy, restartsize_copy);

    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
    dmrginp.set_fullrestart() = false;    

//  // DMRG-LRT: turned off
//  SweepDmrgLrt::do_one(sweepParams, false, direction, false, 0, dmrginp.nroots());
//  sweep_copy.savestate(direction_copy, restartsize_copy);

    break;
  }

  return 0;

}
void calldmrg_(char* input, char* output) {
   int a;
   //a=calldmrg("dmrg.inp",0);//, output);
   a=calldmrg(input,0);//, output);
}

void fullrestartGenblock() {
  SweepParams sweepParams;
  bool direction; int restartsize;
  sweepParams.restorestate(direction, restartsize);
  sweepParams.set_sweep_iter() = 0;
  restartsize = 0;

  SweepGenblock::do_one(sweepParams, false, !direction, RESTART, restartsize, 0);
  
  sweepParams.restorestate(direction, restartsize);
  sweepParams.set_sweep_iter()=0;
  sweepParams.set_block_iter() = 0;
  
  sweepParams.savestate(direction, restartsize);
}  

void restart(double sweep_tol, bool reset_iter)
{
  double last_fe = 100.;
  double last_be = 100.;
  double old_fe = 0.;
  double old_be = 0.;
  int ls_count=0;
  bool direction;
  int restartsize;
  SweepParams sweepParams;
  sweepParams.ls_dw.resize(0);
  sweepParams.ls_energy.resize(0);
  int old_states=sweepParams.get_keep_states();
  int new_states;
  double old_error=0.0;
  double old_energy=0.0;
  bool dodiis = false;

  int domoreIter = 2;

  sweepParams.restorestate(direction, restartsize);
  if(reset_iter) { //this is when you restart from the start of the sweep
    sweepParams.set_sweep_iter() = 0;
    sweepParams.set_restart_iter() = 0;
  }
  
  if (restartwarm)
    last_fe = Sweep::do_one(sweepParams, true, direction, true, restartsize);
  else
    last_fe = Sweep::do_one(sweepParams, false, direction, true, restartsize);


  while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol) || 
	 (dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter()+1 >= sweepParams.get_sweep_iter()) )
  {

    old_fe = last_fe;
    old_be = last_be;
    if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
      break;
    last_be = Sweep::do_one(sweepParams, false, !direction, false, 0);


    //For obtaining the extrapolated energy
    old_states=sweepParams.get_keep_states();
    new_states=sweepParams.get_keep_states_ls();
    if (old_states != new_states) 
    {
       old_energy = last_be+dmrginp.get_coreenergy();
       old_error = sweepParams.get_largest_dw();
       sweepParams.ls_dw.push_back(old_error);
       sweepParams.ls_energy.push_back(old_energy);
       ls_count++;

       if (ls_count >=3) {
          least_squares(sweepParams.ls_dw, sweepParams.ls_energy);
       }
    }


    if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
      break;
    last_fe = Sweep::do_one(sweepParams, false, direction, false, 0);


  }

  if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()){
    //For obtaining the extrapolated energy
    old_energy = last_be+dmrginp.get_coreenergy();
    old_error = sweepParams.get_largest_dw();
    sweepParams.ls_dw.push_back(old_error);
    sweepParams.ls_energy.push_back(old_energy);
    ls_count++;
    if (ls_count >=3) {
       least_squares(sweepParams.ls_dw, sweepParams.ls_energy);
    }
#ifndef MOLPRO
    pout << "Maximum sweep iterations achieved " << std::endl;
#else
    xout << "Maximum sweep iterations achieved " << std::endl;
#endif
  }

}

void dmrg(double sweep_tol)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  int ls_count=0;
  SweepParams sweepParams;
  sweepParams.ls_dw.resize(0);
  sweepParams.ls_energy.resize(0);
  int old_states=sweepParams.get_keep_states();
  int new_states;
  double old_error=0.0;
  double old_energy=0.0;
  // warm up sweep ...
  bool dodiis = false;

  int domoreIter = 0;

  //initialize array of size m_maxiter or dmrginp.max_iter() for dw and energy


  last_fe = Sweep::do_one(sweepParams, true, true, false, 0);
  while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol) || 
	 (dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter()+1 >= sweepParams.get_sweep_iter()) )
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = Sweep::do_one(sweepParams, false, false, false, 0);
      if (dmrginp.outputlevel() > 0) 
         pout << "Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;

      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
         break;

      //For obtaining the extrapolated energy
      old_states=sweepParams.get_keep_states();
      new_states=sweepParams.get_keep_states_ls();
      if (old_states != new_states) 
      {
         old_energy = last_be+dmrginp.get_coreenergy();
         old_error = sweepParams.get_largest_dw();
         sweepParams.ls_dw.push_back(old_error);
         sweepParams.ls_energy.push_back(old_energy);
         ls_count++;

         if (ls_count >=3) {
            least_squares(sweepParams.ls_dw, sweepParams.ls_energy);
         }
      }

      last_fe = Sweep::do_one(sweepParams, false, true, false, 0);

      new_states=sweepParams.get_keep_states();


      if (dmrginp.outputlevel() > 0)
         pout << "Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      if (domoreIter == 2) {
	dodiis = true;
	break;
      }

    }
    if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()) {

    //For obtaining the extrapolated energy
    old_energy = last_be+dmrginp.get_coreenergy();
    old_error = sweepParams.get_largest_dw();
    sweepParams.ls_dw.push_back(old_error);
    sweepParams.ls_energy.push_back(old_energy);
    ls_count++;
    if (ls_count >=3) {
       least_squares(sweepParams.ls_dw, sweepParams.ls_energy);
    }

    pout << "Maximum sweep iterations achieved " << std::endl;
  }

  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());
}

