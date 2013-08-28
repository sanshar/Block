///*                                                                           
//Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
//Copyright (c) 2012, Garnet K.-L. Chan                                        
//                                                                             
//This program is integrated in Molpro with the permission of 
//Sandeep Sharma and Garnet K.-L. Chan
//*/
//
//
//#include "npdm_block_and_decimate.h"
//#include "sweeptwopdm.h"
//#include "global.h"
//#include "solver.h"
//#include "initblocks.h"
//#include "rotationmat.h"
//#include "davidson.h"
//#include "linear.h"
//#include "guess_wavefunction.h"
//#include "twopdm_driver.h"
////#include "threepdm.h"
//#include "density.h"
//#include "davidson.h"
//#include "pario.h"
//
//#ifndef SERIAL
//#include <boost/mpi/communicator.hpp>
//#include <boost/mpi.hpp>
//#endif
//using namespace boost;
//using namespace std;
//
//namespace SpinAdapted {
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//void calcenergy(array_4d<double>& twopdm, int state)
//{
//  using namespace SpinAdapted;
//  Matrix onepdm(2*dmrginp.last_site(), 2*dmrginp.last_site()); onepdm = 0.0;
//  int nelec = dmrginp.real_particle_number();
//
//  for (int i=0; i<dmrginp.last_site()*2; i++)
//  for (int j=0; j<dmrginp.last_site()*2; j++)
//  for (int k=0; k<dmrginp.last_site()*2; k++)
//    onepdm(i+1, j+1) += twopdm(i, k, k, j);
//
//  onepdm /= (nelec-1);
//  double nel = 0.0, sz=0.0;
//  for (int i=0; i<dmrginp.last_site(); i++) {
//    nel += onepdm(2*i+1, 2*i+1)+onepdm(2*i+2, 2*i+2);
//    sz += onepdm(2*i+1, 2*i+1)-onepdm(2*i+2, 2*i+2);
//  }
//
//  double energy = 0.0;
//  for (int i=0; i<dmrginp.last_site()*2; i++)
//  for (int j=0; j<dmrginp.last_site()*2; j++)
//  for (int k=0; k<dmrginp.last_site()*2; k++)
//  for (int l=0; l<dmrginp.last_site()*2; l++)
//    energy += v_2(i,j,k,l)*twopdm(i,j,l,k);
//
//  energy *= 0.5;
//
//  for (int i=0; i<dmrginp.last_site()*2; i++)
//  for (int j=0; j<dmrginp.last_site()*2; j++)
//    energy += v_1(i,j) * onepdm(i+1,j+1);
//
//  pout << "energy of state "<< state <<" = "<< energy+dmrginp.get_coreenergy()<<endl;
//
//  ofstream out("onepdm_fromtpdm");
//  for (int i=0; i<dmrginp.last_site()*2; i++)
//  for (int j=0; j<dmrginp.last_site()*2; j++)
//    out<<i<<"  "<<j<<"  "<<onepdm(i+1,j+1)<<endl;
//  out.close();
//  
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//double SweepTwopdm::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int state)
//{
//  cout.precision(12);
//  SpinBlock system;
//  const int nroots = dmrginp.nroots();
//  std::vector<double> finalEnergy(nroots,0.);
//  std::vector<double> finalEnergy_spins(nroots,0.);
//  double finalError = 0.;
//
//  sweepParams.set_sweep_parameters();
//  // a new renormalisation sweep routine
//  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
//  pout << "\t\t\t ============================================================================ " << endl;
//  
//  InitBlocks::InitStartingBlock (system,forward, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp);
//  if(!restart)
//    sweepParams.set_block_iter() = 0;
// 
//  pout << "\t\t\t Starting block is :: " << endl << system << endl;
//  if (!restart) 
//    SpinBlock::store (forward, system.get_sites(), system); // if restart, just restoring an existing block --
//  sweepParams.savestate(forward, system.get_sites().size());
//  bool dot_with_sys = true;
//
//  //MAW
//  Twopdm_driver twopdm_driver;
//  twopdm_driver.resize_array(2*dmrginp.last_site());
//  twopdm_driver.clear_array();
//  for (int i=0; i<nroots; i++)
//    twopdm_driver.save_npdm_binary(i, i); 
//
//
//  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
//    {
//      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
//      pout << "\t\t\t ----------------------------" << endl;
//      if (forward)
//	pout << "\t\t\t Current direction is :: Forwards " << endl;
//      else
//	pout << "\t\t\t Current direction is :: Backwards " << endl;
//
//      //if (SHOW_MORE) pout << "system block" << endl << system << endl;
//  
//      if (dmrginp.no_transform())
//	      sweepParams.set_guesstype() = BASIC;
//      else if (!warmUp && sweepParams.get_block_iter() != 0) 
//  	    sweepParams.set_guesstype() = TRANSFORM;
//      else if (!warmUp && sweepParams.get_block_iter() == 0 && 
//                ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
//                  dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
//        sweepParams.set_guesstype() = TRANSPOSE;
//      else
//        sweepParams.set_guesstype() = BASIC;
//      
//      pout << "\t\t\t Blocking and Decimating " << endl;
//	  
//      SpinBlock newSystem;
//
//      // Build 2pdm elements
//      Npdm::BlockAndDecimate(twopdm_driver, sweepParams, system, newSystem, warmUp, dot_with_sys, state);
//
//      for(int j=0;j<nroots;++j)
//        pout << "\t\t\t Total block energy for State [ " << j << 
//	  " ] with " << sweepParams.get_keep_states()<<" :: " << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              
//
//      finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
//      finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
//      finalError = max(sweepParams.get_lowest_error(),finalError);
//
//      system = newSystem;
//
//      pout << system<<endl;
//      
//      SpinBlock::store (forward, system.get_sites(), system);	 	
//
//      pout << "\t\t\t saving state " << system.get_sites().size() << endl;
//      ++sweepParams.set_block_iter();
//      sweepParams.savestate(forward, system.get_sites().size());
//    }
//  //for(int j=0;j<nroots;++j)
//  {int j = state;
//    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
//	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
//  }
//  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
//  pout << "\t\t\t ============================================================================ " << endl;
//
//
//  int i = state, j = state;
//  twopdm_driver.load_npdm_binary(i, j); 
////MAW >>>>>>>>>>>>>>>
////print twopdm
////std::cout << "Final 2PDM:\n";
////for (int i=0; i<2*dmrginp.last_site(); ++i) {
////for (int j=0; j<2*dmrginp.last_site(); ++j) {
////for (int k=0; k<2*dmrginp.last_site(); ++k) {
////for (int l=0; l<2*dmrginp.last_site(); ++l) {
////  if ( abs(twopdm(i,j,k,l)) > 1e-12 ) 
////    std::cout << "maw-so-2pdm  " << i << "," << j << "," << k << "," << l << "\t\t" << twopdm(i,j,k,l) << std::endl;
////}   
////}
////}
////}
////MAW <<<<<<<<<<<<<
//  //calcenergy(twopdm, i);
//  twopdm_driver.save_npdm_text(i, j);
//  twopdm_driver.save_spatial_npdm_text(i, j);
//  twopdm_driver.save_spatial_npdm_binary(i, j);
//  
//
//  // update the static number of iterations
//
//  ++sweepParams.set_sweep_iter();
//
//  return finalEnergy[0];
//
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//}
