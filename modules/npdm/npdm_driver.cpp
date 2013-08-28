/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <fstream>
// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

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

unsigned int get_mpi_tag( int rank0, int rank1, int lda )
{
  return 100 * (rank0 * lda + rank1);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::pair<bool, NpdmSpinOps_base> >
Npdm_driver::get_all_mpi_ops( const bool local_skip, NpdmSpinOps & local_ops, std::vector< boost::mpi::request > & reqs )
{
  boost::mpi::communicator world;
  std::vector< std::pair<bool, NpdmSpinOps_base> > all_ops;
  reqs.clear();

  // First element is local set of spin operators
  NpdmSpinOps_base local_base(local_ops);
  std::pair<bool, NpdmSpinOps_base> local_pair = std::make_pair( local_skip, local_base );
  all_ops.push_back( local_pair );

  // Serial calculation
  if (world.size() == 1) return all_ops;

//--- boost serialization broken for NpdmSpinOps_base object!!!   
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//        // save data to archive
//    boost::shared_ptr<SparseMatrix> op (new Cre);
//    op = local_base.opReps_.at(0);
//          std::ofstream ofs("crap.tmp");
//          boost::archive::text_oarchive oa(ofs);
//          oa << *(local_base.opReps_.at(0));
//          oa << *op;
//          ofs.close();
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


  // Do MPI blocking communication //FIXME non-blocking?
  for (int rank = 0; rank < world.size(); ++rank) {
    if ( rank != mpigetrank() ) {
      // Send to non-local rank
      unsigned int tag = get_mpi_tag(mpigetrank(), rank, world.size());
      assert( tag%100 == 0 );
      world.send(rank, tag, local_skip);
      local_base.send_mpi_obj(rank, tag+1, tag+50);
      // Recv from non-local rank
      bool nonlocal_skip;
      NpdmSpinOps_base nonlocal_base;
      tag = get_mpi_tag(rank, mpigetrank(), world.size());
      world.recv(rank, tag, nonlocal_skip);
      nonlocal_base.recv_mpi_obj(rank, tag+1, tag+50, local_base.opReps_.size());
      // Store non-local data
      std::pair<const bool, NpdmSpinOps_base> nonlocal_pair = std::make_pair( nonlocal_skip, nonlocal_base );
      all_ops.push_back( nonlocal_pair );
    }
  }

//  // Nonlocal is other in 2-proc case
//  bool nonlocal_skip;
//  NpdmSpinOps_base nonlocal_base;
//
//  if ( mpigetrank() == 0) {
//    // Send to 1
//    world.send(1, 0, local_skip);
//    local_base.send_mpi_obj(1, 200, 250);
//    // Recv from 1
//    world.recv(1, 100, nonlocal_skip);
//    nonlocal_base.recv_mpi_obj(1, 300, 350, local_base.opReps_.size());
//  }
//  else if ( mpigetrank() == 1) {
//    // Send to 0
//    world.send(0, 100, local_skip);
//    local_base.send_mpi_obj(0, 300, 350);
//    // Recv from 0
//    world.recv(0, 0, nonlocal_skip);
//    nonlocal_base.recv_mpi_obj(0, 200, 250, local_base.opReps_.size());
//  }
//  else
//    assert(false);
//
//  std::pair<const bool, NpdmSpinOps_base> nonlocal_pair = std::make_pair( nonlocal_skip, nonlocal_base );
//  all_ops.push_back( nonlocal_pair );

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

  // All MPI threads must be synchronised here so they all work on same operator pattern simultaneously
  std::cout.flush();
  world.barrier();
  int lhs_maxsize = get_mpi_max_lhs_size( lhsOps.size() );

  // Only one spatial combination on the dot block (including NULL)
  assert( dotOps.size() == 1 );
  bool skip = dotOps.set_local_ops( 0 );
//FIXME is this skip OK in parallel?
  if (skip) return;

  // Many spatial combinations on left block
  assert( lhsOps.size() > 0 );
  for ( int ilhs = 0; ilhs < lhs_maxsize; ++ilhs ) {
    // Set local operators as dummy if load-balancing isn't perfect
    if ( ilhs < lhsOps.size() )
      skip = lhsOps.set_local_ops( ilhs );
    else
      skip = true;

    // Collect LHS ops across all MPI ranks
    std::vector< boost::mpi::request > reqs;
    std::vector< std::pair<bool, NpdmSpinOps_base> > all_lhsOps = get_all_mpi_ops( skip, lhsOps, reqs ); 

//cout << "about to call loop over LHS; skip = " << skip << " ; rank = " << mpigetrank() <<  std::endl;
    // Contract all LHS ops with local RHS ops
    for ( auto lhs_mpi_ops = all_lhsOps.begin(); lhs_mpi_ops != all_lhsOps.end(); ++lhs_mpi_ops ) {
      if ( ! lhs_mpi_ops->first ) {
//cout << "about to call inner loop over LHS; rank = " << mpigetrank() <<  std::endl;
        do_npdm_inner_loop( npdm_expectations, lhs_mpi_ops->second, rhsOps, dotOps ); 
      }
    }
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

//===========================================================================================================================================================

}


  //    if (skip) continue;
  
  //    std::vector< boost::shared_ptr<NpdmSpinOps> > all_lhsOps;
  //    boost::mpi::all_gather( world, lhsOps, all_lhsOps );
  //    assert( all_lhsOps.size() == world.size() );
  
  
  //cout << "DEBUG doing!\n";
  //DEBUG >>>>>>>>>>>>>>>>>>>>>>
  //BROKEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //////    // save data to archive
  //////        std::ofstream ofs("crap.tmp");
  //////        boost::archive::text_oarchive oa(ofs);
  //////
  //////        oa.register_type<Npdm_op_wrapper_compound_CCDD>();
  //////
  //////
  //////        // write class instance to archive
  //////cout << "MAW about to serialize\n";
  //////        oa << lhsOps;
  //////        ofs.close();
  //////
  //////cout << "MAW about to read back\n";
  //////    // ... some time later restore the class instance to its orginal state
  ////////    boost::shared_ptr<NpdmSpinOps> op_new (new Npdm_op_wrapper_compound_CCDD(lhsBlock));
  ////////    boost::shared_ptr<NpdmSpinOps> op_new;
  //////      Npdm_op_wrapper_compound_CCDD op_new(lhsBlock);
  //////cout << "MAW 1\n";
  //////        // create and open an archive for input
  //////        std::ifstream ifs("crap.tmp");
  //////cout << "MAW 2\n";
  //////        boost::archive::text_iarchive ia(ifs);
  //////cout << "MAW 3\n";
  //////        // read class state from archive
  //////        ia >> op_new;
  //////cout << "op_new.factor_ , lhsOps->factor_  = " << op_new.factor_ << "," << lhsOps->factor_ << std::endl;
  //////assert(op_new.factor_ == lhsOps->factor_);
  //////assert(op_new.transpose_ == lhsOps->transpose_);
  //////        // archive and stream closed when destructors are called
  //////        
  ///--------------------------------------------------------------------------------------
  
  ////  ONLY NEED TO COMMUNICATE SKIP, INDICES_ AND OPREPS_, EVERYTHING ELSE SHOULD BE THE SAME!
  //        boost::shared_ptr<NpdmSpinOps> op_new (new Npdm_op_wrapper_compound_CCDD(lhsBlock));
  //    // save data to archive
  //        std::ofstream ofs("crap.tmp");
  //        boost::archive::text_oarchive oa(ofs);
  //        oa << lhsOps->indices_;
  //        ofs.close();
  //        // communciate
  //        std::ifstream ifs("crap.tmp");
  //        boost::archive::text_iarchive ia(ifs);
  //        ia >> op_new->indices_;
  //      assert( op_new->indices_.size() == lhsOps->indices_.size() );
  //      cout << "op_new->indices_[0] , lhsOps->indices_[0]  = " << op_new->indices_[0] << "," << lhsOps->indices_[0] << std::endl;
  //      assert(op_new->indices_ == lhsOps->indices_);
  //
  ///// OPREPS
  //        assert ( lhsOps->opReps_.size() > 0 );
  //        // create and open an archive for opreps
  //        std::ofstream ofs2("crap2.tmp");
  //        boost::archive::text_oarchive oa2(ofs2);
  //   cout << "about to write opReps! size = " << lhsOps->opReps_.size() << std::endl;
  //        oa2 << *(lhsOps->opReps_.at(0));
  //        ofs2.close();
  //   cout << "done writing opReps!\n";
  //        // communciate
  //        std::ifstream ifs2("crap2.tmp");
  //        boost::archive::text_iarchive ia2(ifs2);
  //   cout << "about to read opReps!  size = " << op_new->opReps_.size() << std::endl;
  //    boost::shared_ptr<SparseMatrix> op (new Cre);
  ////    ia2 >> *(lhsOps->opReps_.at(0));
  //    ia2 >> *op;
  //    assert( op_new->opReps_.size() == 0 );
  //    op_new->opReps_.push_back(op) ;
  //    cout << "lhsOps:\n";
  //    cout << *lhsOps->opReps_.at(0) << std::endl;
  //    cout << "new_op:\n";
  //    cout << *op_new->opReps_.at(0) << std::endl;
  //
  //
  //cout << "DEBUG tested!\n";
  ////assert(false);
  //DEBUG <<<<<<<<<<<<<<<<<<<<<<
  
  
//      std::vector< boost::mpi::request > reqs;
//      std::vector< std::pair<bool, boost::shared_ptr<NpdmSpinOps>> > all_lhsOps = get_all_mpi_ops( skip, lhsOps, reqs ); 
  
//      for ( auto lhs_mpi_op = all_lhsOps.begin(); lhs_mpi_op != all_lhsOps.end(); ++lhs_mpi_op ) {
//        // First element is local lhsOps so no need to wait
//        if ( ! lhs_mpi_op->first ) {
//          do_npdm_inner_loop( npdm_expectations, lhs_mpi_op->second, rhsOps, dotOps ); 
//        }
//        // Check all MPI communication of lhsOps is finished (or part of it?) and wait if necessary
//  //std::cout << "Waiting for MPI comm to finish... " << mpigetrank() << std::endl;
//  //      boost::mpi::wait_all(reqs.begin(), reqs.end());
//  //std::cout << "MPI comm done!  " << mpigetrank() << std::endl;
//  //      reqs.clear();
//      }
  
//    // Synchronize all MPI ranks here ??
//    }
//  
//    // Close file if needed (put in wrapper destructor??)
//    if ( lhsOps.ifs_.is_open() ) lhsOps.ifs_.close();
//  
//  }
//  

//----------------------------------------------
//
//  boost::shared_ptr<NpdmSpinOps> nonlocal_ops (new NpdmSpinOps);
////FIXME what kind of copy/assignment is this???
/////  nonlocal_ops = local_ops;
//  nonlocal_ops->size_ = local_ops->size_;
//  nonlocal_ops->mults_ = local_ops->mults_;
//  nonlocal_ops->build_pattern_ = local_ops->build_pattern_;
//  nonlocal_ops->transpose_ = local_ops->transpose_;
//  nonlocal_ops->factor_ = local_ops->factor_;
//
////  assert( nonlocal_ops->indices_ == local_ops->indices_ );
////  assert( nonlocal_ops->opReps_.size() > 0 );
//// Seems to be a shallow copy!
////  nonlocal_ops->opReps_.clear();
////  assert( local_ops->opReps_.size() > 0 );
////  nonlocal_ops->indices_.clear();
////  assert( local_ops->indices_.size() > 0 );
//
//
////    boost::shared_ptr<SparseMatrix> op (new Cre);
//////    ia2 >> *(lhsOps->opReps_.at(0));
////    ia2 >> *op;
////    assert( op_new->opReps_.size() == 0 );
////    op_new->opReps_.push_back(op) ;
//
//  nonlocal_ops->opReps_.clear();
//  for ( int i = 0; i < local_ops->opReps_.size(); ++i) {
//    boost::shared_ptr<SparseMatrix> op (new Cre);
//    nonlocal_ops->opReps_.push_back(op);
//  }
//
//  // Now send operators to other MPI ranks in non-blocking fashion
//  bool nonlocal_skip;
////cout << "MAW mpi ops " << i << " communicate rank=" << mpigetrank() << std::endl;
//
//  if ( mpigetrank() == 0) {
//    // Send
//    world.send(1, 0, local_skip);
//    world.send(1, 1, local_ops->indices_);
//    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
//      world.send(1, i+2, *(local_ops->opReps_.at(i)) );
////    reqs.push_back( world.isend(1, 0, local_skip) );
////    reqs.push_back( world.isend(1, 1, local_ops->indices_) );
////    reqs.push_back( world.isend(1, 2, local_ops->opReps_) );
//    // Recv
//    world.recv(1, 30, nonlocal_skip);
//    world.recv(1, 40, nonlocal_ops->indices_);
//    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
//      world.recv(1, i+50, *(nonlocal_ops->opReps_.at(i)) );
////    reqs.push_back( world.irecv(1, 3, nonlocal_skip) );
////    reqs.push_back( world.irecv(1, 4, nonlocal_ops->indices_) );
////    reqs.push_back( world.irecv(1, 5, nonlocal_ops->opReps_) );
//  }
//  else if ( mpigetrank() == 1) {
//    // Send
//    world.send(0, 30, local_skip);
//    world.send(0, 40, local_ops->indices_);
//    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
//      world.send(0, i+50, *(local_ops->opReps_.at(i)) );
////    reqs.push_back( world.isend(0, 3, local_skip) );
////    reqs.push_back( world.isend(0, 4, local_ops->indices_) );
////    reqs.push_back( world.isend(0, 5, local_ops->opReps_) );
//    // Recv
//    world.recv(0, 0, nonlocal_skip);
//    world.recv(0, 1, nonlocal_ops->indices_);
//    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
//      world.recv(0, i+2, *(nonlocal_ops->opReps_.at(i)) );
////    reqs.push_back( world.irecv(0, 0, nonlocal_skip) );
////    reqs.push_back( world.irecv(0, 1, nonlocal_ops->indices_) );
////    reqs.push_back( world.irecv(0, 2, nonlocal_ops->opReps_) );
//  }
//  else
//    assert(false);
//
//cout << local_ops->build_pattern_ << "; local skip = " << local_skip << std::endl;
//cout << nonlocal_ops->build_pattern_ << "; nonlocal skip = " << nonlocal_skip << std::endl;
//  std::pair<bool, boost::shared_ptr<NpdmSpinOps>> nonlocal_pair = std::make_pair( nonlocal_skip, nonlocal_ops );
//  all_ops.push_back( nonlocal_pair );
//
//
//
//  return all_ops;
//}
//  
