/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm_block_and_decimate.h"
#include "sweepfourpdm.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include "fourpdm.h"
#include "density.h"
#include "davidson.h"
#include "pario.h"

#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
using namespace boost;
using namespace std;


namespace SpinAdapted{

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double SweepFourpdm::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int state)
{
  cout.precision(12);
  SpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
  if(!restart)
    sweepParams.set_block_iter() = 0;
 
  pout << "\t\t\t Starting block is :: " << endl << system << endl;
  if (!restart) 
    SpinBlock::store (forward, system.get_sites(), system); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

  int dim = 2*dmrginp.last_site();
  array_8d<double> fourpdm(dim,dim,dim,dim,dim,dim,dim,dim);
  fourpdm.Clear();
  for (int i=0; i<nroots; i++)
    save_fourpdm_binary(fourpdm, i, i); 


  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	pout << "\t\t\t Current direction is :: Forwards " << endl;
      else
	pout << "\t\t\t Current direction is :: Backwards " << endl;

      //if (SHOW_MORE) pout << "system block" << endl << system << endl;
  
      if (dmrginp.no_transform())
	      sweepParams.set_guesstype() = BASIC;
      else if (!warmUp && sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if (!warmUp && sweepParams.get_block_iter() == 0 && 
                ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
                  dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;
      
      pout << "\t\t\t Blocking and Decimating " << endl;
	  
      SpinBlock newSystem;

      Npdm::BlockAndDecimate ("fourpdm", sweepParams, system, newSystem, warmUp, dot_with_sys, state);

      for(int j=0;j<nroots;++j)
        pout << "\t\t\t Total block energy for State [ " << j << 
	  " ] with " << sweepParams.get_keep_states()<<" :: " << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              

      finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
      finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
      finalError = max(sweepParams.get_lowest_error(),finalError);

      system = newSystem;

      pout << system<<endl;
      
      SpinBlock::store (forward, system.get_sites(), system);	 	

      pout << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      sweepParams.savestate(forward, system.get_sites().size());
    }
  //for(int j=0;j<nroots;++j)
  {int j = state;
    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
  }
  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t ============================================================================ " << endl;


  int i = state, j = state;
  //for (int j=0; j<=i; j++) {
  load_fourpdm_binary(fourpdm, i, j); 
std::cout << "4PDM done!\n";
std::cout << "WARNING: not saving full 4PDM spatial binary!\n";
  save_fourpdm_text(fourpdm, i, j);
  save_spatial_fourpdm_text(fourpdm, i, j);
//  save_spatial_fourpdm_binary(fourpdm, i, j);
  

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalEnergy[0];
//}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
