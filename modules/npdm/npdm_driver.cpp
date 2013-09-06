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

namespace SpinAdapted{

// DEBUG only
std::vector<double> DEBUG_COMM_TIME(1000);
std::vector<int> DEBUG_CALL_GET_EXPECT(1000);

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================
//void check_operator_count( NpdmSpinOps& Ops )
//{
//  boost::mpi::communicator world;
//
//  // Serial case
//  if ( world.size() == 1 ) return;
//
//  pout << "Finding duplicate ops\n";
//  std::vector< int > indices = Ops.get_1d_indices();
//
//  if ( mpigetrank() == 0 ) {
//    // Gather indices from all ranks
//    std::vector< std::vector<int> > all_indices;
//    boost::mpi::gather( world, indices, all_indices, 0 );
//    // Identify duplicates
//    std::set<int> unique;
//    std::set<int> duplicates;
//    for (int rank=0; rank < world.size(); ++rank ) { 
//      for (int k=0; k < all_indices.at(rank).size(); ++k ) { 
//        int id = all_indices[rank][k];
//        if (unique.find(id) != unique.end() )
//          unique.insert(id);
//        else
//          duplicates.insert(id);
//      }
//    }
//    cout << "Duplicates:\n";
//    for (auto it = duplicates.begin(); it != duplicates.end(); ++it) {
//      cout << *it << endl;
//    }
//  }
//  else {
//    gather( world, indices, 0 );
//  }
//}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::save_full_array(int i, int j) 
{
  Timer timer;
  // Combine NPDM elements from all mpi ranks and dump to file
  accumulate_npdm();
  save_spatial_npdm_text(i, j);
//FIXME  save_npdm_binary(i, j);
//FIXME  save_npdm_text(i, j);
//FIXME  save_spatial_npdm_binary(i, j);

  pout << "NPDM save full array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::save_sparse_array(int i, int j) 
{
boost::mpi::communicator world;
world.barrier();
  Timer timer;
  // Save nonredundant npdm elements from each mpi rank as separate files
  sparse_array_.dump_file(i,j);
//FIXME
world.barrier();
  pout << "NPDM save sparse array time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

unsigned int get_mpi_tag( int rank0, int rank1, int lda )
{
  unsigned int tag = rank0 * lda + rank1;
  assert( tag < 42949672 );
  return 100 * tag;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

int Npdm_driver::get_mpi_max_size( int my_size )
{
  int maxsize;
  boost::mpi::communicator world;
  std::vector<int> all_sizes;
  boost::mpi::all_gather(world, my_size, all_sizes);
  maxsize = *std::max_element( all_sizes.begin(), all_sizes.end() );
  assert( my_size <= maxsize );
  return maxsize;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Do we parallelize by broadcasting LHS or RHS operators?  This is a very simple heuristic for now.  With disk-access, it may not be so good.

bool Npdm_driver::broadcast_lhs( int lhs_size, int rhs_size )
{
  // Note all ranks have to make the same decision!
  bool do_lhs = true;
  int lhs_maxsize = get_mpi_max_size( lhs_size );
  int rhs_maxsize = get_mpi_max_size( rhs_size );
  if (rhs_maxsize < lhs_maxsize) do_lhs = false;
  return do_lhs;
}
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::do_inner_loop( const char inner, Npdm::Npdm_expectations& npdm_expectations, 
                                 NpdmSpinOps_base& outerOps, NpdmSpinOps& innerOps, NpdmSpinOps& dotOps ) 
{

  if ( dotOps.opReps_.size() > 0 ) assert( dotOps.mults_.size() == dotOps.opReps_.size() );
  if ( outerOps.opReps_.size() > 0 ) assert( outerOps.mults_.size() == outerOps.opReps_.size() );

  // Many spatial combinations on right block
  for ( int iop = 0; iop < innerOps.size(); ++iop ) {
    bool skip = innerOps.set_local_ops( iop );
    if (skip) continue;
    if ( innerOps.opReps_.size() > 0 ) assert( innerOps.mults_.size() == innerOps.opReps_.size() );

    // Get non-spin-adapated 3PDM elements after building spin-adapted elements
    std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements;
    DEBUG_CALL_GET_EXPECT[mpigetrank()] += 1;
    // This always works out as calling in order (lhs,rhs,dot); or at least it should!
    if ( inner == 'r' )
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( outerOps, innerOps, dotOps );
    else if ( inner == 'l' )
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( innerOps, outerOps, dotOps );
    else
      assert(false);

    // Store new npdm elements
    sparse_array_.insert( new_spin_orbital_elements );
    if ( use_full_array_ ) assign_npdm_elements( new_spin_orbital_elements );
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Originally wrote this routine for LHS and RHS operators, but also works interchanging "lhs" and "rhs" everywhere when called in reverse

void Npdm_driver::loop_over_block_operators( const char inner, Npdm::Npdm_expectations & npdm_expectations,
                                             NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps, bool lhsdot ) 
{
  boost::mpi::communicator world;
  int lhs_maxsize = get_mpi_max_size( lhsOps.size() );

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhs_maxsize; ++ilhs ) {
    // Set local operators as dummy if load-balancing isn't perfect
    bool skip = true;
    if ( ilhs < lhsOps.size() ) skip = lhsOps.set_local_ops( ilhs );

    // Don't want to parallelize if any of these conditions hold (either makes no sense or creates duplicated work)
    bool serial = ( world.size() == 1            ||
                    lhsOps.build_pattern_ == ""  ||
                    rhsOps.build_pattern_ == ""  ||
                    lhsdot );

    if ( serial ) {
      if ( !skip ) do_inner_loop( inner, npdm_expectations, lhsOps, rhsOps, dotOps ); 
    }
    else {
      // Parallelize by broadcasting LHS ops
      NpdmSpinOps_base local_base(lhsOps);
      std::vector< boost::mpi::request > reqs;
      std::vector< NpdmSpinOps_base > nonlocal_base( world.size() );

      Timer timer;
      // Communicate array sizes etc (note blocking until comm is complete)
      std::vector< int > nonlocal_size( world.size() );
      std::vector< int > nonlocal_skip( world.size() );
      int local_skip = skip; // apparently boost::mpi::all_gather fails with bools...??
      boost::mpi::all_gather(world, local_skip, nonlocal_skip);
      int local_size = local_base.opReps_.size();
      boost::mpi::all_gather(world, local_size, nonlocal_size);
      DEBUG_COMM_TIME[mpigetrank()] += timer.elapsedwalltime();
    
      // Communicate operator reps
      for (int rank = 0; rank < world.size(); ++rank) {
        if ( rank != mpigetrank() ) {
          // Get unique tag for send-recv pair (asymmetric)
          unsigned int send_tag = get_mpi_tag(mpigetrank(), rank, world.size()); assert( send_tag%100 == 0 );
          unsigned int recv_tag = get_mpi_tag(rank, mpigetrank(), world.size()); assert( send_tag%100 == 0 );
          if ( ! local_skip ) {
            std::vector< boost::mpi::request > new_reqs = local_base.isend_mpi_obj(rank, send_tag+2, send_tag+50);
            reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );
          }
          if ( ! nonlocal_skip.at(rank) ) {
            std::vector< boost::mpi::request > new_reqs = nonlocal_base.at(rank).irecv_mpi_obj(rank, recv_tag+2, recv_tag+50, nonlocal_size.at(rank));
            reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );
          }
        }
      }
      
      // Do loop over RHS with local LHS operator while waiting for all non-local to be communicated (FIXME can we extend this idea?)
      if ( ! local_skip ) do_inner_loop( inner, npdm_expectations, local_base, rhsOps, dotOps ); 

      // Contract all nonlocal LHS ops with local RHS ops; must wait for communication to be finished first
      Timer timer2;
      boost::mpi::wait_all( reqs.begin(), reqs.end() );
      DEBUG_COMM_TIME[mpigetrank()] += timer2.elapsedwalltime();
      for (int rank = 0; rank < world.size(); ++rank) {
        if ( rank != mpigetrank() ) {
          if ( ! nonlocal_skip.at(rank) ) do_inner_loop( inner, npdm_expectations, nonlocal_base.at(rank), rhsOps, dotOps ); 
        }
      }

      // Synchronize all MPI ranks here
      Timer timer3;
      std::cout.flush();
      world.barrier();
      DEBUG_COMM_TIME[mpigetrank()] += timer3.elapsedwalltime();
    }
  }

  // Close file if needed (put in wrapper destructor??)
  if ( lhsOps.ifs_.is_open() ) lhsOps.ifs_.close();

}
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos)
{
  boost::mpi::communicator world;
  std::cout.flush();
  world.barrier();
  DEBUG_COMM_TIME[mpigetrank()] = 0;
  DEBUG_CALL_GET_EXPECT[mpigetrank()] = 0;
  Timer timer;
  pout << "===========================================================================================\n";
  pout << "NPDM sweep position = "<< sweepPos << " (" << endPos << ")\n";

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

    // MPI threads must be synchronised here so they all work on same operator pattern simultaneously
    Timer timer3;
    std::cout.flush();
    world.barrier();
    DEBUG_COMM_TIME[mpigetrank()] += timer3.elapsedwalltime();

    //pout << "-------------------------------------------------------------------------------------------\n";
    //pout << "Doing pattern:  ";
    npdm_patterns.print_cd_string( pattern->at('l') );
    npdm_patterns.print_cd_string( pattern->at('d') );
    npdm_patterns.print_cd_string( pattern->at('r') );
    pout << std::endl;
    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');

    // Choice of read from disk or not done inside the wrapper
    boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
    boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
    boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );
    //cout << "lhsOps.size() = " << lhsOps->size() << "; rank = " << mpigetrank() <<  std::endl;
    //cout << "rhsOps.size() = " << rhsOps->size() << "; rank = " << mpigetrank() <<  std::endl;

    // Only one spatial combination on the dot block (including NULL)
    assert( dotOps->size() == 1 );
    bool skip = dotOps->set_local_ops( 0 );
    if ( ! skip ) {
      // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
      if ( broadcast_lhs( lhsOps->size(), rhsOps->size() ) ) {
        loop_over_block_operators( 'r', npdm_expectations, *lhsOps, *rhsOps, *dotOps, (lhsBlock->size() == 1) );
      }
      else {
        loop_over_block_operators( 'l', npdm_expectations, *rhsOps, *lhsOps, *dotOps, (rhsBlock->size() == 1) );
      }
    }
  }
  if (world.rank() == 0) {
    double sum;
    reduce(world, DEBUG_COMM_TIME[mpigetrank()], sum, std::plus<double>(), 0);
    pout << "NPDM mpi communications time " << sum << endl;
  } else {
    reduce(world, DEBUG_COMM_TIME[mpigetrank()], std::plus<double>(), 0);
  }
  if (world.rank() == 0) {
    int sum;
    reduce(world, DEBUG_CALL_GET_EXPECT[mpigetrank()], sum, std::plus<int>(), 0);
    pout << "NPDM calls to expectation engine " << sum << endl;
  } else {
    reduce(world, DEBUG_CALL_GET_EXPECT[mpigetrank()], std::plus<int>(), 0);
  }
  pout << "NPDM compute elements time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
  
}

//===========================================================================================================================================================

}
