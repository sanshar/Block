/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <boost/mpi.hpp>
#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "npdm_driver.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"

#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "rotationmat.h"
#include "davidson.h"
#include "linear.h"
#include "guess_wavefunction.h"
#include "density.h"
#include "davidson.h"
#include "pario.h"

namespace SpinAdapted{

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================

unsigned int get_mpi_tag( int rank0, int rank1, int lda )
{
  unsigned int tag = rank0 * lda + rank1;
  assert( tag < 42949672 );
  return 100 * tag;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector<NpdmSpinOps_base> Npdm_driver::get_all_mpi_ops( const bool local_skip, NpdmSpinOps & local_ops, std::vector< boost::mpi::request > & reqs )
{
  boost::mpi::communicator world;
  std::vector< NpdmSpinOps_base > all_ops;
  reqs.clear();

  // First element is local set of spin operators
  NpdmSpinOps_base local_base(local_ops);
  if ( ! local_skip ) all_ops.push_back( local_base );

  // Serial calculation
  if (world.size() == 1) return all_ops;

  // Do MPI blocking communication //FIXME non-blocking?
  for (int rank = 0; rank < world.size(); ++rank) {
    if ( rank != mpigetrank() ) {

      // Send to non-local rank
      unsigned int tag = get_mpi_tag(mpigetrank(), rank, world.size());
      assert( tag%100 == 0 );
      world.send(rank, tag, local_skip);
      int local_size = local_base.opReps_.size();
      world.send(rank, tag+1, local_size);
      if ( ! local_skip ) local_base.send_mpi_obj(rank, tag+2, tag+50);

      // Recv from non-local rank
      bool nonlocal_skip;
      NpdmSpinOps_base nonlocal_base;
      tag = get_mpi_tag(rank, mpigetrank(), world.size());
      world.recv(rank, tag, nonlocal_skip);
      int nonlocal_size;
      world.recv(rank, tag+1, nonlocal_size);
      if ( ! nonlocal_skip ) nonlocal_base.recv_mpi_obj(rank, tag+2, tag+50, nonlocal_size);

      // Store non-local data
      if ( ! nonlocal_skip ) all_ops.push_back( nonlocal_base );
    }
  }

  return all_ops;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::do_npdm_inner_loop( Npdm::Npdm_expectations & npdm_expectations, NpdmSpinOps_base & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps ) 
{

  if ( dotOps.opReps_.size() > 0 ) assert( dotOps.mults_.size() == dotOps.opReps_.size() );
  if ( lhsOps.opReps_.size() > 0 ) assert( lhsOps.mults_.size() == lhsOps.opReps_.size() );

  // Many spatial combinations on right block
  for ( int irhs = 0; irhs < rhsOps.size(); ++irhs ) {
    bool skip = rhsOps.set_local_ops( irhs );
    if (skip) continue;
    if ( rhsOps.opReps_.size() > 0 ) assert( rhsOps.mults_.size() == rhsOps.opReps_.size() );

    // Get non-spin-adapated 3PDM elements after building spin-adapted elements
    std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements;
    new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( lhsOps, rhsOps, dotOps );

    // Assign npdm elements
    assign_npdm_elements( new_spin_orbital_elements );
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
int Npdm_driver::get_mpi_max_lhs_size( int my_size )
{
  int maxsize;
  boost::mpi::communicator world;
  std::vector<int> all_sizes;
  all_gather(world, my_size, all_sizes);
  maxsize = *std::max_element( all_sizes.begin(), all_sizes.end() );
  assert( my_size <= maxsize );
  return maxsize;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::npdm_loop_over_block_operators( Npdm::Npdm_expectations & npdm_expectations,
                                                  NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps ) 
{

  boost::mpi::communicator world;
//cout << "-------------------------------------------------------------------------------------------\n";
cout << "lhsOps.size() = " << lhsOps.size() << "; rank = " << mpigetrank() <<  std::endl;
//cout << "dotOps.size() = " << dotOps.size() << "; rank = " << mpigetrank() <<  std::endl;
//cout << "rhsOps.size() = " << rhsOps.size() << "; rank = " << mpigetrank() <<  std::endl;

  // MPI threads must be synchronised here so they all work on same operator pattern simultaneously
std::cout.flush();
  world.barrier();
  int lhs_maxsize = get_mpi_max_lhs_size( lhsOps.size() );

  // Only one spatial combination on the dot block (including NULL)
  assert( dotOps.size() == 1 );
  bool skip = dotOps.set_local_ops( 0 );
//FIXME is this skip OK in parallel?
  if (skip) return;

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhs_maxsize; ++ilhs ) {
    // Set local operators as dummy if load-balancing isn't perfect
    if ( ilhs < lhsOps.size() )
      skip = lhsOps.set_local_ops( ilhs );
    else
      skip = true;

    // Collect LHS ops across all MPI ranks
    std::vector< boost::mpi::request > reqs;
    std::vector< NpdmSpinOps_base > all_lhsOps = get_all_mpi_ops( skip, lhsOps, reqs ); 

    // Contract all LHS ops with local RHS ops
    for ( auto lhs_mpi_ops = all_lhsOps.begin(); lhs_mpi_ops != all_lhsOps.end(); ++lhs_mpi_ops ) 
      do_npdm_inner_loop( npdm_expectations, *lhs_mpi_ops, rhsOps, dotOps ); 

    // Synchronize all MPI ranks here
    std::cout.flush();
    world.barrier();
  }

  // Close file if needed (put in wrapper destructor??)
  if ( lhsOps.ifs_.is_open() ) lhsOps.ifs_.close();

}
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::compute_npdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos)
{

  boost::mpi::communicator world;
  pout << "===========================================================================================\n";
  pout << "NPDM sweep position = "<< sweepPos << " (" << endPos << ")\n";

  // Npdm array built so far
  load_npdm_binary(state, state);
  
  // Initialize class that computes expectation values when sent LHS, Dot and RHS operator spin-sets from this spin-block
  Npdm::Npdm_expectations npdm_expectations( npdm_order_, wavefunctions.at(0), big );

  // Get LHS, Dot and RHS spin-blocks
  SpinBlock* rhsBlock = big.get_rightBlock();
  SpinBlock* lhsdotBlock = big.get_leftBlock();
  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();

  // Loop over NPDM operator patterns
  Npdm::Npdm_patterns npdm_patterns( npdm_order_, sweepPos, endPos );

  for (auto pattern = npdm_patterns.ldr_cd_begin(); pattern != npdm_patterns.ldr_cd_end(); ++pattern) {
//    pout << "===========================================================================================\n";
    cout << "=============================================================================== " << mpigetrank() << std::endl;

    std::cout << "Doing pattern: rank " << mpigetrank() << std::endl;
    npdm_patterns.print_cd_string( pattern->at('l') );
    npdm_patterns.print_cd_string( pattern->at('d') );
    npdm_patterns.print_cd_string( pattern->at('r') );
    std::cout << std::endl;
    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');

    // CHOICE OF READ FROM DISK OR NOT DONE INSIDE THE WRAPPER!!
    boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
    boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
    boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );

//FIXME rewind files if wrapper reads from disk on RHS
//    assert( rhs_cd_type.size() < 3 );

//FIXME
//if ( lhs_cd_type.size() > 3 ) continue;
//if ( lhs_cd_type.size() < 2 ) continue;
//FIXME
//std::vector<Npdm::CD> lhs = {Npdm::CREATION,Npdm::DESTRUCTION,Npdm::CREATION};
//if ( (lhs_cd_type != lhs) ) continue;
//std::vector<Npdm::CD> rhs = {Npdm::CREATION,Npdm::DESTRUCTION,Npdm::DESTRUCTION};
//if ( (rhs_cd_type != rhs) ) continue;
//if ( (lhs_cd_type == foo) || (rhs_cd_type == foo) || (dot_cd_type == foo) ) continue;
//if ( lhs_cd_type.size() == 2  &&
//     dot_cd_type.size() == 2  &&
//     rhs_cd_type.size() == 2 )
    // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
    npdm_loop_over_block_operators( npdm_expectations, *lhsOps, *rhsOps, *dotOps );
  }
  
  // Combine NPDM elements from this sweep point with others
  accumulate_npdm();
  save_npdm_binary(state, state);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::npdm_block_and_decimate( SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, 
                                           const bool &useSlater, const bool& dot_with_sys, int state)
{
  //mcheck("at the start of block and decimate");
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  SpinBlock systemDot;
  SpinBlock envDot;
  int systemDotStart, systemDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  if (forward)
  {
    systemDotStart = *system.get_sites().rbegin () + 1;
    systemDotEnd = systemDotStart + systemDotSize;
  }
  else
  {
    systemDotStart = system.get_sites() [0] - 1;
    systemDotEnd = systemDotStart - systemDotSize;
  }
  vector<int> spindotsites(2); 
  spindotsites[0] = systemDotStart;
  spindotsites[1] = systemDotEnd;
  //if (useSlater) {
    systemDot = SpinBlock(systemDotStart, systemDotEnd);
    //SpinBlock::store(true, systemDot.get_sites(), systemDot);
    //}
    //else
    //SpinBlock::restore(true, spindotsites, systemDot);
  SpinBlock environment, environmentDot, newEnvironment;

  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addAdditionalCompOps();
//FIXME MAW change depending on forward or backward which operators are assigned to which mpi procs
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.get_sys_add(), dmrginp.direct(), DISTRIBUTED_STORAGE, true, true);
  
  InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, true, true, true);
  SpinBlock big;
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 

  const int nroots = dmrginp.nroots();
  std::vector<Wavefunction> solution(1);

  DiagonalMatrix e;
  GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, true, 0.0); 

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, solution, 0);
#endif

  std::vector<Matrix> rotateMatrix;
  DensityMatrix tracedMatrix;
  tracedMatrix.allocate(newSystem.get_stateInfo());
  tracedMatrix.makedensitymatrix(solution, big, std::vector<double>(1,1.0), 0.0, 0.0, false);
  rotateMatrix.clear();
  if (!mpigetrank())
    double error = newSystem.makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  

#ifndef SERIAL
  mpi::broadcast(world,rotateMatrix,0);
#endif

//MAW
  int sweepPos = sweepParams.get_block_iter();
  int endPos = sweepParams.get_n_iters()-1;
  compute_npdm_sweep(solution, big, state, sweepPos, endPos);
//MAW

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);

  //for(int i=0;i<dmrginp.nroots();++i)
  solution[0].SaveWavefunctionInfo (big.get_stateInfo(), big.get_leftBlock()->get_sites(), state);

  newSystem.transform_operators(rotateMatrix);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double Npdm_driver::do_one_sweep(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int state)
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

  //MAW
  npdm_resize_array(2*dmrginp.last_site());
  npdm_clear_array();
  for (int i=0; i<nroots; i++)
    save_npdm_binary(i, i); 


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

      // Build Npdm elements
      npdm_block_and_decimate(sweepParams, system, newSystem, warmUp, dot_with_sys, state);

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
  load_npdm_binary(i, j); 
//MAW >>>>>>>>>>>>>>>
//print twopdm
//std::cout << "Final 2PDM:\n";
//for (int i=0; i<2*dmrginp.last_site(); ++i) {
//for (int j=0; j<2*dmrginp.last_site(); ++j) {
//for (int k=0; k<2*dmrginp.last_site(); ++k) {
//for (int l=0; l<2*dmrginp.last_site(); ++l) {
//  if ( abs(twopdm(i,j,k,l)) > 1e-12 ) 
//    std::cout << "maw-so-2pdm  " << i << "," << j << "," << k << "," << l << "\t\t" << twopdm(i,j,k,l) << std::endl;
//}   
//}
//}
//}
//MAW <<<<<<<<<<<<<
  //calcenergy(twopdm, i);
  save_npdm_text(i, j);
  save_spatial_npdm_text(i, j);
  save_spatial_npdm_binary(i, j);
  

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalEnergy[0];

}

//===========================================================================================================================================================

}
