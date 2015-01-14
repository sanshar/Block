/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "global.h"
#include "Symmetry.h"
#include "sweep_params.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

using namespace boost;

SpinAdapted::SweepParams::SweepParams()
{
  //calculationType = DMRG;
  restart_iter = 0;
  block_iter = 0;
  sweep_iter = 0;
  n_iters = 0;
  error = 0.0;
  largest_dw = 0.0;
  forward_starting_size = 1;
  backward_starting_size = 1;
  keep_states = 0;
  keep_qstates = 0;
  noise = 0.0;
  additional_noise = 0.0;
  davidson_tol = 1.e-6;
  guesstype = BASIC;

  onedot = (dmrginp.algorithm_method() == ONEDOT);
  sys_add = dmrginp.sys_add();
  env_add = dmrginp.env_add();
  forward_starting_size = 1;
  if(dmrginp.spinAdapted())
    n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
  else
    n_iters = (dmrginp.last_site() - 4*forward_starting_size - 2*sys_add - 2*env_add) / (2*sys_add) + 1;

  backward_starting_size = forward_starting_size;
}

void SpinAdapted::SweepParams::set_sweep_parameters()
{
  int iter;
  int current = 0;
  //ROA: Deactivated least squares
  //int current_ls = 0;
  int sweep_iter_ls = sweep_iter + 1;
  for (iter = 0; iter < dmrginp.sweep_iter_schedule().size(); ++iter)
  {
     if (sweep_iter >= dmrginp.sweep_iter_schedule()[iter]) { 
       current = iter; 
       //pout << "chosen current " << iter << endl;
     }
  }
  /*
  //ROA: Deactivated least squares
  for (iter = 0; iter < dmrginp.sweep_iter_schedule().size(); ++iter)
  {
     if (sweep_iter_ls >= dmrginp.sweep_iter_schedule()[iter]) { 
       current_ls = iter; 
       //pout << "chosen current " << iter << endl;
     }
  }
  */
  keep_states = dmrginp.sweep_state_schedule()[current];
  //ROA: Deactivated least squares
  //keep_states_ls = dmrginp.sweep_state_schedule()[current_ls];
  keep_qstates = 0.0;//dmrginp.sweep_qstate_schedule()[current];
  davidson_tol = dmrginp.sweep_tol_schedule()[current];
           
  if (dmrginp.get_twodot_method() == 1) {
   if (this->get_sweep_iter() == 0) 
    this->set_additional_noise() = dmrginp.get_twodot_noise();
   if (this->get_sweep_iter() > 0 && (this->get_sweep_iter() < dmrginp.twodot_to_onedot_iter() || dmrginp.algorithm_method()==TWODOT))
    this->set_additional_noise() = 0.5*dmrginp.get_twodot_gamma()*this->get_largest_dw();
  }

  noise = dmrginp.sweep_noise_schedule()[current];
  additional_noise = this->get_additional_noise();//dmrginp.get_twodot_noise();

   p1out << "\n\t\t\t Sweep iteration ... " << SpinAdapted::SweepParams::sweep_iter;
   p1out<< " "  << SpinAdapted::SweepParams::keep_states << " " << SpinAdapted::SweepParams::davidson_tol << " " << SpinAdapted::SweepParams::noise << " "<< SpinAdapted::SweepParams::additional_noise <<endl; 
  else 
    pout << endl;

  if (dmrginp.algorithm_method() == ONEDOT) {
    onedot = true;
    env_add = 0;
    if(dmrginp.spinAdapted())
      n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
    else
      n_iters = (dmrginp.last_site() - 4*forward_starting_size - 2*sys_add - 2*env_add) / (2*sys_add) + 1;
    //n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
    pout << "\t\t\t Using the one dot algorithm ... " << endl;
  }  
  if(dmrginp.algorithm_method() == TWODOT_TO_ONEDOT) {
    if(dmrginp.twodot_to_onedot_iter() <= SpinAdapted::SweepParams::sweep_iter)
    {
      onedot = true;
      env_add = 0;
      if(dmrginp.spinAdapted())
	n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
      else
	n_iters = (dmrginp.last_site() - 4*forward_starting_size - 2*sys_add - 2*env_add) / (2*sys_add) + 1;
      if (dmrginp.twodot_to_onedot_iter() == SpinAdapted::SweepParams::sweep_iter)
	pout << "\t\t\t Switching from two dot to one dot ... " << endl;
    }
    else {
      onedot = false;
      env_add = 1;
      if(dmrginp.spinAdapted())
	n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
      else
	n_iters = (dmrginp.last_site() - 4*forward_starting_size - 2*sys_add - 2*env_add) / (2*sys_add) + 1;
    }
  }
  if(dmrginp.algorithm_method() == TWODOT) {
    onedot = false;
    env_add = 1;
    if(dmrginp.spinAdapted())
      n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
    else
      n_iters = (dmrginp.last_site() - 4*forward_starting_size - 2*sys_add - 2*env_add) / (2*sys_add) + 1;
	pout << "\t\t\t Using the two dot algorithm ... " << endl;
  }

  //now figure out number of iterations and starting size, during first call only
  p2out << "\t\t\t forward system starting size ... " << forward_starting_size << " " << n_iters << endl;
  p2out << "\t\t\t backward system starting size ... " << backward_starting_size << " " << n_iters << endl;
  //pout << "onedot or twodot "<<dmrginp.algorithm_method()<<endl;

}


void SpinAdapted::SweepParams::calc_niter()
{
    if(dmrginp.spinAdapted())
      n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
    else
      n_iters = (dmrginp.last_site() - 4*forward_starting_size - 2*sys_add - 2*env_add) / (2*sys_add) + 1;
    //n_iters = (dmrginp.last_site() - 2*forward_starting_size - sys_add - env_add) / sys_add + 1;
}

void SpinAdapted::SweepParams::savestate(const bool &forward, const int &size)
{
  if (mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(), "/statefile.", mpigetrank(), ".tmp");
    p1out << "\t\t\t Saving state "<<file<<endl;
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
    p1out << "\t\t\t Loading state "<<file<<endl;
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load_wave(ifs);
    load_wave >> forward >> size >> *this;
    ifs.close();
  }
  size = block_iter+1;
  restart_iter = sweep_iter;

  if (block_iter >= get_n_iters()) {
    sweep_iter ++;
    block_iter = 0;
    size = block_iter +1;
    forward = !forward;
  }


  pout << "\n\t\t\t Restarting at sweep iteration: " << sweep_iter <<endl;
  pout << "\n\t\t\t Restarting at block iterations: "  << block_iter <<endl;

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, forward, 0); 
  mpi::broadcast(world, size, 0);
  mpi::broadcast(world, *this, 0);
#endif
}

