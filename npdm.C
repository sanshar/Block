/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm.h"
#include "sweepgenblock.h"
#include "sweep_params.h"
#include "sweeponepdm.h"
#include "twopdm_driver.h"
#include "threepdm_driver.h"
#include "fourpdm_driver.h"

void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
void ReadInput(char* conf);
void fullrestartGenblock();

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm( int npdm_order )
{

  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;

  if (dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm not allowed with twodot algorithm" << endl;
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
  dmrginp.do_npdm_ops() = true;
  dmrginp.screen_tol() = 0.0;

  sweep_copy.restorestate(direction_copy, restartsize_copy);
  dmrginp.set_fullrestart() = true;
  sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
//FIXME don't compute unnecessary operators until NOW
  SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
  dmrginp.set_fullrestart() = false;

  switch (npdm_order) {
  case (1):
    // Compute onepdm elements
    SweepOnepdm::do_one(sweepParams, false, direction, false, 0);
    sweep_copy.savestate(direction_copy, restartsize_copy);
    break;
  case (2):
    // Compute twopdm elements
    for (int state=0; state<dmrginp.nroots(); state++) {
      Twopdm_driver twopdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      twopdm_driver.do_one_sweep(sweepParams, false, direction, false, 0, state);
    }
    break;
  case (3):
    // Compute threepdm elements
    for (int state=0; state<dmrginp.nroots(); state++) {
      Threepdm_driver threepdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      threepdm_driver.do_one_sweep(sweepParams, false, direction, false, 0, state);
    }
    break;
  case (4):
    // Compute fourpdm elements
    for (int state=0; state<dmrginp.nroots(); state++) {
      Fourpdm_driver fourpdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      fourpdm_driver.do_one_sweep(sweepParams, false, direction, false, 0, state);
    }
    break;
  }
  sweep_copy.savestate(direction_copy, restartsize_copy);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm_restart( int npdm_order )
{
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;

  if(sym == "dinfh") {
    pout << "Npdm not implemented with dinfh symmetry"<<endl;
    abort();
  }

  sweepParams.restorestate(direction, restartsize);
  if(!sweepParams.get_onedot() || dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm only runs for the onedot algorithm" << endl;
    abort();
  }

  dmrginp.screen_tol() = 0.0; //need to turn screening off for onepdm
  dmrginp.Sz() = dmrginp.total_spin_number();
  dmrginp.do_npdm_ops() = true;
  dmrginp.screen_tol() = 0.0;

  sweep_copy.restorestate(direction_copy, restartsize_copy);
  dmrginp.set_fullrestart() = true;
  sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
  SweepGenblock::do_one(sweepParams, false, !direction, false, 0, 0); //this will generate the cd operators
  dmrginp.set_fullrestart() = false;

  switch (npdm_order) {
  case (1):
    SweepOnepdm::do_one(sweepParams, false, direction, false, 0);
    sweep_copy.savestate(direction_copy, restartsize_copy);
    break;
  case (2):
    for (int state=0; state<dmrginp.nroots(); state++) {
      Twopdm_driver twopdm_driver;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      twopdm_driver.do_one_sweep(sweepParams, false, direction, false, 0, state);
    }
    break;
  }
  sweep_copy.savestate(direction_copy, restartsize_copy);

}

//===========================================================================================================================================================

}
}

