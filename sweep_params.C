#include "global.h"
#include "Symmetry.h"
#include "sweep_params.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace boost;

SpinAdapted::SweepParams::SweepParams()
{
  //calculationType = DMRG;
  restart_iter = 0;
  block_iter = 0;
  sweep_iter = 0;
  n_iters = 0;
  error = 0.0;
  forward_starting_size = 1;
  backward_starting_size = 1;
  keep_states = 0;
  keep_qstates = 0;
  noise = 0;
  additional_noise = 0;
  davidson_tol = 1.e-6;
  guesstype = BASIC;

  onedot = (dmrginp.algorithm_method() == ONEDOT);
  sys_add = dmrginp.sys_add();
  env_add = dmrginp.env_add();
  forward_starting_size = 1;
  n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
  backward_starting_size = forward_starting_size;
}

void SpinAdapted::SweepParams::set_sweep_parameters()
{
  int iter;
  int current = 0;
  for (iter = 0; iter < dmrginp.sweep_iter_schedule().size(); ++iter)
  {
     if (sweep_iter >= dmrginp.sweep_iter_schedule()[iter]) { 
       current = iter; 
       //pout << "chosen current " << iter << endl;
     }
  }
  keep_states = dmrginp.sweep_state_schedule()[current];
  keep_qstates = 0.0;//dmrginp.sweep_qstate_schedule()[current];
  davidson_tol = dmrginp.sweep_tol_schedule()[current];
  noise = dmrginp.sweep_noise_schedule()[current];
  additional_noise = 0.0;//dmrginp.sweep_additional_noise_schedule()[current];

  if (dmrginp.outputlevel() != 0) {
   pout << "\t\t\t Sweep iteration ... " << SpinAdapted::SweepParams::sweep_iter;
   pout<< " "  << SpinAdapted::SweepParams::keep_states << " " << SpinAdapted::SweepParams::davidson_tol << " " << SpinAdapted::SweepParams::noise << " "<< SpinAdapted::SweepParams::additional_noise <<endl; 
  }
  else 
    pout << endl;
  //now figure out number of iterations and starting size, during first call only
  if (dmrginp.outputlevel() != 0) {
    pout << "\t\t\t forward system starting size ... " << forward_starting_size << " " << n_iters << endl;
    pout << "\t\t\t backward system starting size ... " << backward_starting_size << " " << n_iters << endl;
    //pout << "onedot or twodot "<<dmrginp.algorithm_method()<<endl;
  }
  if (dmrginp.algorithm_method() == ONEDOT) {
    onedot = true;
    env_add = 0;
    //n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
    pout << "\t\t\t Using the one dot algorithm ... " << endl;
  }  
  if(dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() <= SpinAdapted::SweepParams::sweep_iter)
  {
    onedot = true;
    env_add = 0;
    n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
    if (dmrginp.twodot_to_onedot_iter() == SpinAdapted::SweepParams::sweep_iter)
      pout << "\t\t\t Switching from two dot to one dot ... " << endl;
  }
}


void SpinAdapted::SweepParams::calc_niter()
{
  n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
}

void SpinAdapted::SweepParams::savestate(const bool &forward, const int &size)
{
  if (mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(), "/statefile.", mpigetrank(), ".tmp");
    if (dmrginp.outputlevel() != 0)
      pout << "\t\t\t Saving state "<<file<<endl;
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save_wave(ofs);
    save_wave << forward << size << *this;
    ofs.close();
  }
}

void SpinAdapted::SweepParams::restorestate(bool &forward, int &size)
{
  if (mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d%s", dmrginp.load_prefix().c_str(), "/statefile.", mpigetrank(), ".tmp");
    if (dmrginp.outputlevel() != 0)
      pout << "\t\t\t Loading state "<<file<<endl;
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load_wave(ifs);
    load_wave >> forward >> size >> *this;
    ifs.close();
  }
  size = block_iter+1;
  restart_iter = sweep_iter;

  /*
  if (set_calcType() == FCI || set_calcType() == TINYCALC) {
    if (dmrginp.calc_type() != DMRG) 
      pout << "Restart not supported after a FCI  or a TINYCALC calculation.";
    else 
      pout << "Cannot perform Genblock/Onepdm/Twopdm calculation after a FCI or a TINYCALC calculation."<<endl;
    abort();
  }


  if (dmrginp.algorithm_method() == TWODOT) 
    onedot = false;
  env_add = dmrginp.env_add();
  forward_starting_size = 1;
  n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
  backward_starting_size = forward_starting_size;
  */
  if (block_iter >= get_n_iters()) {
    sweep_iter ++;
    block_iter = 0;
    size = block_iter +1;
    forward = !forward;
  }


  pout << "\t\t\t Restarting at sweep iteration: " << sweep_iter <<endl;
  pout << "\t\t\t Restarting at block iterations: "  << block_iter <<endl;

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, forward, 0); 
  mpi::broadcast(world, size, 0);
  mpi::broadcast(world, *this, 0);
#endif
}

