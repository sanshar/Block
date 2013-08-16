/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <boost/format.hpp>
//FIXME#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/export.hpp>
//#include <boost/serialization/utility.hpp> // for std::pair
//DEBUG >>>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
//DEBUG <<<
//FIXME#endif
#include "execinfo.h"

#include "twopdm.h"
#include "npdm_driver.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"
#include "npdm_epermute.h"

namespace SpinAdapted{

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================
// FIXME put in header instead?
Npdm_driver::Npdm_driver(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big) 
: wavefunctions_(wavefunctions), big_(big)
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::do_npdm_inner_loop( array_4d<double> & npdm, 
                                      boost::shared_ptr<NpdmSpinOps> & lhsOps,
                                      boost::shared_ptr<NpdmSpinOps> & rhsOps,
                                      boost::shared_ptr<NpdmSpinOps> & dotOps,
                                      Npdm::Npdm_expectations & npdm_expectations )
{

//pout << "lhsOps->indices_.size() " << lhsOps->indices_.size() << std::endl;
//pout << "lhsOps->opReps_.size() " << lhsOps->opReps_.size() << std::endl;
//pout << "dotOps->indices_.size() " << dotOps->indices_.size() << std::endl;
//pout << "dotOps->opReps_.size() " << dotOps->opReps_.size() << std::endl;
    if ( lhsOps->opReps_.size() > 0 ) assert( lhsOps->mults_.size() == lhsOps->opReps_.size() );

    // Many spatial combinations on right block
    for ( int irhs = 0; irhs < rhsOps->size(); ++irhs ) {
      bool skip = rhsOps->set_local_ops( irhs );
      if (skip) continue;
//pout << "rhsOps->indices_.size() " << rhsOps->indices_.size() << std::endl;
//pout << "rhsOps->opReps_.size() " << rhsOps->opReps_.size() << std::endl;
//pout << "-------------------------------------------------------------------------------------------\n";
//pout << "spatial: ilhs, irhs = " << ilhs << ", " << irhs << std::endl;
      if ( rhsOps->opReps_.size() > 0 ) assert( rhsOps->mults_.size() == rhsOps->opReps_.size() );

      // Get non-spin-adapated 3PDM elements after building spin-adapted elements
      std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( 6 );

      // Assign npdm elements
      assign_twopdm_elements( new_spin_orbital_elements, npdm );
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::pair<bool, boost::shared_ptr<NpdmSpinOps>> > 
Npdm_driver::get_all_mpi_ops( const bool local_skip, boost::shared_ptr<NpdmSpinOps>& local_ops, std::vector< boost::mpi::request > & reqs )
{
  boost::mpi::communicator world;
  std::vector< std::pair<bool, boost::shared_ptr<NpdmSpinOps>> > all_ops;
  reqs.clear();

  // First element is local set of spin operators
  std::pair<bool, boost::shared_ptr<NpdmSpinOps>> local_pair = std::make_pair( local_skip, local_ops );
  all_ops.push_back( local_pair );
//FIXME only distribute ops that are non-local
//if ( local_ops->indices_.size() > 3 ) assert(false);
//if ( local_ops->indices_.size() < 2 ) assert(false);
return all_ops;
//----------------------------------------------

  boost::shared_ptr<NpdmSpinOps> nonlocal_ops (new NpdmSpinOps);
//FIXME what kind of copy/assignment is this???
///  nonlocal_ops = local_ops;
  nonlocal_ops->size_ = local_ops->size_;
  nonlocal_ops->mults_ = local_ops->mults_;
  nonlocal_ops->build_pattern_ = local_ops->build_pattern_;
  nonlocal_ops->transpose_ = local_ops->transpose_;
  nonlocal_ops->factor_ = local_ops->factor_;

//  assert( nonlocal_ops->indices_ == local_ops->indices_ );
//  assert( nonlocal_ops->opReps_.size() > 0 );
// Seems to be a shallow copy!
//  nonlocal_ops->opReps_.clear();
//  assert( local_ops->opReps_.size() > 0 );
//  nonlocal_ops->indices_.clear();
//  assert( local_ops->indices_.size() > 0 );


//    boost::shared_ptr<SparseMatrix> op (new Cre);
////    ia2 >> *(lhsOps->opReps_.at(0));
//    ia2 >> *op;
//    assert( op_new->opReps_.size() == 0 );
//    op_new->opReps_.push_back(op) ;

  nonlocal_ops->opReps_.clear();
  for ( int i = 0; i < local_ops->opReps_.size(); ++i) {
    boost::shared_ptr<SparseMatrix> op (new Cre);
    nonlocal_ops->opReps_.push_back(op);
  }

  // Now send operators to other MPI ranks in non-blocking fashion
  bool nonlocal_skip;
//  for (int i=0; i < world.size(); ++i) {
//cout << "MAW mpi ops " << i << " communicate rank=" << mpigetrank() << std::endl;

  if ( mpigetrank() == 0) {
    // Send
    world.send(1, 0, local_skip);
    world.send(1, 1, local_ops->indices_);
    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
      world.send(1, i+2, *(local_ops->opReps_.at(i)) );
//    reqs.push_back( world.isend(1, 0, local_skip) );
//    reqs.push_back( world.isend(1, 1, local_ops->indices_) );
//    reqs.push_back( world.isend(1, 2, local_ops->opReps_) );
    // Recv
    world.recv(1, 30, nonlocal_skip);
    world.recv(1, 40, nonlocal_ops->indices_);
    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
      world.recv(1, i+50, *(nonlocal_ops->opReps_.at(i)) );
//    reqs.push_back( world.irecv(1, 3, nonlocal_skip) );
//    reqs.push_back( world.irecv(1, 4, nonlocal_ops->indices_) );
//    reqs.push_back( world.irecv(1, 5, nonlocal_ops->opReps_) );
  }
  else if ( mpigetrank() == 1) {
    // Send
    world.send(0, 30, local_skip);
    world.send(0, 40, local_ops->indices_);
    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
      world.send(0, i+50, *(local_ops->opReps_.at(i)) );
//    reqs.push_back( world.isend(0, 3, local_skip) );
//    reqs.push_back( world.isend(0, 4, local_ops->indices_) );
//    reqs.push_back( world.isend(0, 5, local_ops->opReps_) );
    // Recv
    world.recv(0, 0, nonlocal_skip);
    world.recv(0, 1, nonlocal_ops->indices_);
    for ( int i = 0; i < local_ops->opReps_.size(); ++i)
      world.recv(0, i+2, *(nonlocal_ops->opReps_.at(i)) );
//    reqs.push_back( world.irecv(0, 0, nonlocal_skip) );
//    reqs.push_back( world.irecv(0, 1, nonlocal_ops->indices_) );
//    reqs.push_back( world.irecv(0, 2, nonlocal_ops->opReps_) );
  }
  else
    assert(false);

cout << local_ops->build_pattern_ << "; local skip = " << local_skip << std::endl;
cout << nonlocal_ops->build_pattern_ << "; nonlocal skip = " << nonlocal_skip << std::endl;
  std::pair<bool, boost::shared_ptr<NpdmSpinOps>> nonlocal_pair = std::make_pair( nonlocal_skip, nonlocal_ops );
  all_ops.push_back( nonlocal_pair );



  return all_ops;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::npdm_loop_over_block_operators( std::vector<Npdm::CD> & lhs_cd_type,
                                                  std::vector<Npdm::CD> & dot_cd_type,
                                                  std::vector<Npdm::CD> & rhs_cd_type,
                                                  array_4d<double> & npdm )
{
  boost::mpi::communicator world;

  SpinBlock* rhsBlock = big_.get_rightBlock();
  SpinBlock* lhsdotBlock = big_.get_leftBlock();
  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();

// CHOICE OF READ FROM DISK OR NOT DONE INSIDE THE WRAPPER!!
  boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
  boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );

  Npdm::Npdm_expectations npdm_expectations( wavefunctions_.at(0), big_, *lhsOps, *dotOps, *rhsOps );
//cout << "-------------------------------------------------------------------------------------------\n";
cout << "lhsOps->size() = " << lhsOps->size() << "; rank = " << mpigetrank() <<  std::endl;
//cout << "dotOps->size()" << dotOps->size() << std::endl;
//cout << "rhsOps->size()" << rhsOps->size() << std::endl;
//cout << "------\n";

  // Only one spatial combination on the dot block
  assert( dotOps->size() == 1 );
  bool skip = dotOps->set_local_ops( 0 );
  if (skip) return;
  if ( lhsOps->opReps_.size() > 0 ) assert( dotOps->mults_.size() == dotOps->opReps_.size() );

  // Many spatial combinations on left block
  assert( lhsOps->size() > 0 );
  for ( int ilhs = 0; ilhs < lhsOps->size(); ++ilhs ) {
    skip = lhsOps->set_local_ops( ilhs );
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


    std::vector< boost::mpi::request > reqs;
    std::vector< std::pair<bool, boost::shared_ptr<NpdmSpinOps>> > all_lhsOps = get_all_mpi_ops( skip, lhsOps, reqs ); 

    for ( auto lhs_mpi_op = all_lhsOps.begin(); lhs_mpi_op != all_lhsOps.end(); ++lhs_mpi_op ) {
      // First element is local lhsOps so no need to wait
      if ( ! lhs_mpi_op->first ) {
        do_npdm_inner_loop( npdm, lhs_mpi_op->second, rhsOps, dotOps, npdm_expectations );
      }
      // Check all MPI communication of lhsOps is finished (or part of it?) and wait if necessary
//std::cout << "Waiting for MPI comm to finish... " << mpigetrank() << std::endl;
//      boost::mpi::wait_all(reqs.begin(), reqs.end());
//std::cout << "MPI comm done!  " << mpigetrank() << std::endl;
//      reqs.clear();
    }

  // Synchronize all MPI ranks here ??
  }

  // Close file if needed (put in wrapper destructor??)
  if ( lhsOps->ifs_.is_open() ) lhsOps->ifs_.close();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::compute_npdm_sweep(int state, int sweepPos, int endPos) {

  pout << "===========================================================================================\n";
  pout << "2PDM sweep position = "<< sweepPos << " (" << endPos << ")\n";

  // 3pdm array built so far
// 2pdm array built so far

  int dim = 2*big_.size();
  array_4d<double> npdm(dim,dim,dim,dim);
  load_twopdm_binary(npdm, state, state);
  
  // Loop over NPDM operator patterns (here we initialize for 3PDM)
// Loop over NPDM operator patterns (here we initialize for 2PDM)
  Npdm::Npdm_patterns npdm_patterns( 2, sweepPos, endPos );

  for (auto pattern = npdm_patterns.ldr_cd_begin(); pattern != npdm_patterns.ldr_cd_end(); ++pattern) {

    pout << "===========================================================================================\n";
    std::cout << "Doing pattern:\n";
    npdm_patterns.print_cd_string( pattern->at('l') );
    npdm_patterns.print_cd_string( pattern->at('d') );
    npdm_patterns.print_cd_string( pattern->at('r') );
    std::cout << std::endl;

    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');

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
    npdm_loop_over_block_operators( lhs_cd_type, dot_cd_type, rhs_cd_type, npdm );
  }
  
  // Combine NPDM elements from this sweep point with others
  accumulate_twopdm(npdm);
  save_twopdm_binary(npdm, state, state);

}

//===========================================================================================================================================================

}
