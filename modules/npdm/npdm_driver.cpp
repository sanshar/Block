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
  save_npdm_binary(i, j);
  // Dump spatial and text NPDM files
  save_npdm_text(i, j);
  save_spatial_npdm_text(i, j);
  save_spatial_npdm_binary(i, j);
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
cout << "dotOps.size() = " << dotOps.size() << "; rank = " << mpigetrank() <<  std::endl;
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
  
}

//===========================================================================================================================================================

}
