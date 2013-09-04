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

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================

void Npdm_driver::save_full_array(int i, int j) 
{
  // Combine NPDM elements from all mpi ranks and save one binary file
  accumulate_npdm();
//FIXME  save_npdm_binary(i, j);
  // Dump spatial and text NPDM files
//FIXME  save_npdm_text(i, j);
  save_spatial_npdm_text(i, j);
//FIXME  save_spatial_npdm_binary(i, j);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::save_sparse_array(int i, int j) 
{
  // Save nonredundant npdm elements from each mpi rank as separate files
  sparse_array_.dump_file(i,j);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

unsigned int get_mpi_tag( int rank0, int rank1, int lda )
{
  unsigned int tag = rank0 * lda + rank1;
  assert( tag < 42949672 );
  return 100 * tag;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//std::vector<NpdmSpinOps_base> Npdm_driver::get_all_mpi_ops( const bool local_skip, NpdmSpinOps& local_ops, std::vector< boost::mpi::request >& reqs )
//{
//  boost::mpi::communicator world;
//  std::vector< NpdmSpinOps_base > all_ops;
//  reqs.clear();
//
//  // First element is local set of spin operators
//  NpdmSpinOps_base local_base(local_ops);
//  if ( ! local_skip ) all_ops.push_back( local_base );
//
//  // Serial calculation
//  if (world.size() == 1) return all_ops;
//
//  // Communicate array sizes etc
//  std::vector< int > nonlocal_size( world.size() );
//  std::vector< int > nonlocal_skip( world.size() );
//  int local_skip_i = local_skip; // apparently boost::mpi::all_gather fails with bools...??
//  boost::mpi::all_gather(world, local_skip_i, nonlocal_skip);
//  int local_size = local_base.opReps_.size();
//  boost::mpi::all_gather(world, local_size, nonlocal_size);
//
//  // Communicate operator reps
//  std::vector< NpdmSpinOps_base > nonlocal_base( world.size() );
//  for (int rank = 0; rank < world.size(); ++rank) {
//    if ( rank != mpigetrank() ) {
//
//      // Get unique tag for send-recv pair (asymmetric)
//      unsigned int send_tag = get_mpi_tag(mpigetrank(), rank, world.size()); assert( send_tag%100 == 0 );
//      unsigned int recv_tag = get_mpi_tag(rank, mpigetrank(), world.size()); assert( send_tag%100 == 0 );
//
//      if ( ! local_skip ) {
//        std::vector< boost::mpi::request > new_reqs = local_base.isend_mpi_obj(rank, send_tag+2, send_tag+50);
//        reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );
//      }
//      if ( ! nonlocal_skip.at(rank) ) {
//        std::vector< boost::mpi::request > new_reqs = nonlocal_base.at(rank).irecv_mpi_obj(rank, recv_tag+2, recv_tag+50, nonlocal_size.at(rank));
//        reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );
//      }
//
//    }
//  }
//
////FIXME can we let local lhs op loop over RHS while waiting for all non-local to be communicated?
//boost::mpi::wait_all( reqs.begin(), reqs.end() );
// 
//  // Store non-local data
//  for (int rank = 0; rank < world.size(); ++rank) {
//    if ( rank != mpigetrank() ) {
//      if ( ! nonlocal_skip.at(rank) ) all_ops.push_back( nonlocal_base.at(rank) );
//    }
//  }
//
//  return all_ops;
//}

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

    // Store new npdm elements
    sparse_array_.insert( new_spin_orbital_elements );
    if ( use_full_array_ ) assign_npdm_elements( new_spin_orbital_elements );
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
//cout << "lhsOps.size() = " << lhsOps.size() << "; rank = " << mpigetrank() <<  std::endl;
//cout << "dotOps.size() = " << dotOps.size() << "; rank = " << mpigetrank() <<  std::endl;
//cout << "rhsOps.size() = " << rhsOps.size() << "; rank = " << mpigetrank() <<  std::endl;

  int lhs_maxsize = get_mpi_max_lhs_size( lhsOps.size() );

  // Only one spatial combination on the dot block (including NULL)
  assert( dotOps.size() == 1 );
  bool skip = dotOps.set_local_ops( 0 );
  if (skip) return;

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhs_maxsize; ++ilhs ) {
    // Set local operators as dummy if load-balancing isn't perfect
    if ( ilhs < lhsOps.size() )
      skip = lhsOps.set_local_ops( ilhs );
    else
      skip = true;

    if (world.size() == 1) {
      // Serial calculation
      if ( !skip ) do_npdm_inner_loop( npdm_expectations, lhsOps, rhsOps, dotOps ); 
    }
    else {
      // Parallelize by broadcasting LHS ops
      NpdmSpinOps_base local_base(lhsOps);
      std::vector< boost::mpi::request > reqs;
      std::vector< NpdmSpinOps_base > nonlocal_base( world.size() );

      // Communicate array sizes etc
      std::vector< int > nonlocal_size( world.size() );
      std::vector< int > nonlocal_skip( world.size() );
      int local_skip = skip; // apparently boost::mpi::all_gather fails with bools...??
      boost::mpi::all_gather(world, local_skip, nonlocal_skip);
      int local_size = local_base.opReps_.size();
      boost::mpi::all_gather(world, local_size, nonlocal_size);
    
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
      
      // Do loop over RHS with local LHS operator while waiting for all non-local to be communicated (FIXME can extend this idea?)
      if ( !local_skip ) do_npdm_inner_loop( npdm_expectations, local_base, rhsOps, dotOps ); 
      // Contract all nonlocal LHS ops with local RHS ops; must wait for communication to be finished
      boost::mpi::wait_all( reqs.begin(), reqs.end() );
      for (int rank = 0; rank < world.size(); ++rank) {
        if ( rank != mpigetrank() ) {
          if ( !nonlocal_skip.at(rank) ) do_npdm_inner_loop( npdm_expectations, nonlocal_base.at(rank), rhsOps, dotOps ); 
        }
      }

      // Synchronize all MPI ranks here
      std::cout.flush();
      world.barrier();
    }
  }

  // Close file if needed (put in wrapper destructor??)
  if ( lhsOps.ifs_.is_open() ) lhsOps.ifs_.close();

}
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos)
{

  boost::mpi::communicator world;
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
    pout << "===========================================================================================\n";

    pout << "Doing pattern:\n";
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

    // MPI threads must be synchronised here so they all work on same operator pattern simultaneously
    std::cout.flush();
    world.barrier();
    // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
    npdm_loop_over_block_operators( npdm_expectations, *lhsOps, *rhsOps, *dotOps );
  }
  
}

//===========================================================================================================================================================

}
