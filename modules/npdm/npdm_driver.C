#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

#include "wavefunction.h"
#include "npdm_driver.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "pario.h"

namespace SpinAdapted{
namespace Npdm{

// DEBUG only
std::vector<double> DEBUG_COMM_TIME(1000);
std::vector<int> DEBUG_CALL_GET_EXPECT(1000);

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( SpinBlock * spinBlock, std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================

#ifndef SERIAL
unsigned int get_mpi_tag( int rank0, int rank1, int lda )
{
  unsigned int tag = rank0 * lda + rank1;
  assert( tag < 42949672 );
  return 100 * tag;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef SERIAL
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
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Do we parallelize by broadcasting LHS or RHS operators?  This is a very simple heuristic for now.  With disk-access, it may not be so good.

#ifndef SERIAL
bool Npdm_driver::broadcast_lhs( int lhs_size, int rhs_size )
{
  // Note all ranks have to make the same decision!
  bool do_lhs = true;
  int lhs_maxsize = get_mpi_max_size( lhs_size );
  int rhs_maxsize = get_mpi_max_size( rhs_size );
  if (rhs_maxsize < lhs_maxsize) do_lhs = false;
  return do_lhs;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef SERIAL
bool Npdm_driver::skip_this_mpi_rank( NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps )
{
  boost::mpi::communicator world;
  bool skip = ( mpigetrank() > 0    && 
                lhsOps.is_local_    && 
                rhsOps.is_local_ );
  return skip;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef SERIAL
bool Npdm_driver::skip_parallel( NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, bool lhsrhsdot )
{
  boost::mpi::communicator world;
  // Don't want to parallelize if any of these conditions hold (either makes no sense or creates duplicated work)
  // Assumes 1-index ops are duplicated on all mpi ranks //FIXME check this?
  // (There might be some overlap in these criteria)
  bool skip = ( world.size() == 1               ||
                lhsOps.is_local_                ||
                rhsOps.is_local_                ||
                lhsrhsdot );
  return skip;
}
#endif
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Originally wrote this routine for LHS and RHS operators, but also works interchanging "lhs" and "rhs" everywhere when called in reverse

#ifndef SERIAL
void Npdm_driver::do_parallel_lhs_loop( const char inner, Npdm::Npdm_expectations & npdm_expectations,
                                        NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps, bool skip )
{
  boost::mpi::communicator world;

  // Parallelize by broadcasting LHS ops
  NpdmSpinOps_base local_base(lhsOps);
  std::vector< boost::mpi::request > reqs;
  std::vector< NpdmSpinOps_base > nonlocal_base( world.size() );

  Timer timer;
  // Communicate basic op info
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
  
  // Do loop over RHS with local LHS operator while waiting for all non-local to be communicated 
  // FIXME Can we extend this idea to do batches while other batches are communicating
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
  pout.flush();
  world.barrier();

}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Originally wrote this routine for LHS and RHS operators, but also works interchanging "lhs" and "rhs" everywhere when called in reverse

#ifndef SERIAL
void Npdm_driver::par_loop_over_block_operators( const char inner, Npdm::Npdm_expectations & npdm_expectations,
                                                 NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps, bool lhsrhsdot ) 
{
  int lhs_maxsize = get_mpi_max_size( lhsOps.size() );

  // Skip parallelization completely if it generates duplicates
  if ( skip_this_mpi_rank( lhsOps, rhsOps ) ) return;

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhs_maxsize; ++ilhs ) {
    // Set local operators as dummy if load-balancing isn't perfect
    bool skip_op = true;
    if ( ilhs < lhsOps.size() ) skip_op = lhsOps.set_local_ops( ilhs );

    if ( skip_parallel( lhsOps, rhsOps, lhsrhsdot ) ) {
      if ( ! skip_op ) do_inner_loop( inner, npdm_expectations, lhsOps, rhsOps, dotOps );
    }
    else {
      // Parallelize by broadcasting LHS ops
      do_parallel_lhs_loop( inner, npdm_expectations, lhsOps, rhsOps, dotOps, skip_op );
    }
  }

  assert( ! lhsOps.ifs_.is_open() );
  assert( ! rhsOps.ifs_.is_open() );
}
#endif
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::do_inner_loop( const char inner, Npdm::Npdm_expectations& npdm_expectations, 
                                 NpdmSpinOps_base& outerOps, NpdmSpinOps& innerOps, NpdmSpinOps& dotOps ) 
{
  // Many spatial combinations on right block
  for ( int iop = 0; iop < innerOps.size(); ++iop ) {
    bool skip = innerOps.set_local_ops( iop );
    if (skip) continue;

    // Get non-spin-adapated spin-orbital 3PDM elements after building spin-adapted elements
    std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements;
    DEBUG_CALL_GET_EXPECT[mpigetrank()] += 1;
    // This should always work out as calling in order (lhs,rhs,dot)
    if ( inner == 'r' )
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( outerOps, innerOps, dotOps );
    else if ( inner == 'l' )
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( innerOps, outerOps, dotOps );
    else
      abort();

    // Store new npdm elements
    if ( new_spin_orbital_elements.size() > 0 ) container_.store_npdm_elements( new_spin_orbital_elements );
  }

  assert( ! innerOps.ifs_.is_open() );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::loop_over_block_operators( Npdm::Npdm_expectations& npdm_expectations, NpdmSpinOps& lhsOps, NpdmSpinOps& rhsOps, NpdmSpinOps& dotOps)
{
  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhsOps.size(); ++ilhs ) {
    // Set local operators as dummy if load-balancing isn't perfect
    bool skip_op = true;
    skip_op = lhsOps.set_local_ops( ilhs );
    if ( ! skip_op ) do_inner_loop( 'r', npdm_expectations, lhsOps, rhsOps, dotOps );
  }

  assert( ! lhsOps.ifs_.is_open() );
  assert( ! rhsOps.ifs_.is_open() );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::loop_over_operator_patterns( Npdm::Npdm_patterns& patterns, Npdm::Npdm_expectations& expectations, const SpinBlock& big )
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  // Get LHS, Dot and RHS spin-blocks
  SpinBlock* rhsBlock = big.get_rightBlock();
  SpinBlock* lhsdotBlock = big.get_leftBlock();
  SpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  SpinBlock* dotBlock = lhsdotBlock->get_rightBlock();

  int count = 0;
  for (auto pattern = patterns.ldr_cd_begin(); pattern != patterns.ldr_cd_end(); ++pattern) {
    count++;
    DEBUG_CALL_GET_EXPECT[mpigetrank()] = 0;

#ifndef SERIAL
    // MPI threads must be synchronised here so they all work on same operator pattern simultaneously
    pout.flush();
    world.barrier();
#endif
    //pout << "-------------------------------------------------------------------------------------------\n";
    //pout << "Doing pattern " << count << " of " << patterns.size() << endl;
    //patterns.print_cd_string( pattern->at('l') );
    //patterns.print_cd_string( pattern->at('d') );
    //patterns.print_cd_string( pattern->at('r') );
    //pout << std::endl; 
    //pout.flush();

    // Choice of read from disk or not done inside the wrapper
    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');


    boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
    boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );
    boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );

    // Only one spatial combination on the dot block (including NULL)
    if(dmrginp.spinAdapted()){
    assert( dotOps->size() == 1 );
    bool skip = dotOps->set_local_ops( 0 );
    if ( ! skip ) {
//pout << "p" << mpigetrank() << ": lhs = " << lhsOps->size() << endl;
//pout << "p" << mpigetrank() << ": rhs = " << rhsOps->size() << endl;
      // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
#ifndef SERIAL
      bool lhs_or_rhs_dot = ( (lhsBlock->size() == 1) || (rhsBlock->size() == 1) );
      if ( broadcast_lhs( lhsOps->size(), rhsOps->size() ) ) {
        par_loop_over_block_operators( 'r', expectations, *lhsOps, *rhsOps, *dotOps, lhs_or_rhs_dot );
      }
      else {
        par_loop_over_block_operators( 'l', expectations, *rhsOps, *lhsOps, *dotOps, lhs_or_rhs_dot );
      }
#else
      loop_over_block_operators( expectations, *lhsOps, *rhsOps, *dotOps );
#endif
    }
    }
    else{
      //bool skip=true; 
      // build all kind of dotOps
      for(int i=0; i< dotOps->size(); i++){
        // if it is valid
        if(!dotOps->set_local_ops( i )){

//pout << "p" << mpigetrank() << ": lhs = " << lhsOps->size() << endl;
//pout << "p" << mpigetrank() << ": rhs = " << rhsOps->size() << endl;
      // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
#ifndef SERIAL
      bool lhs_or_rhs_dot = ( (lhsBlock->size() == 1) || (rhsBlock->size() == 1) );
//pout << "lhs_or_rhs_dot " << lhs_or_rhs_dot << endl;
      if ( broadcast_lhs( lhsOps->size(), rhsOps->size() ) ) {
//pout << "broadcast lhs\n";
//pout.flush();
        par_loop_over_block_operators( 'r', expectations, *lhsOps, *rhsOps, *dotOps, lhs_or_rhs_dot );
      }
      else {
//pout << "broadcast rhs\n";
//pout.flush();
        par_loop_over_block_operators( 'l', expectations, *rhsOps, *lhsOps, *dotOps, lhs_or_rhs_dot );
      }
#else
      loop_over_block_operators( expectations, *lhsOps, *rhsOps, *dotOps );
#endif
      }
      }
    }

  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  pout.flush();
  world.barrier();
#endif
  DEBUG_COMM_TIME[mpigetrank()] = 0;
  DEBUG_CALL_GET_EXPECT[mpigetrank()] = 0;
  Timer timer;
  pout << "===========================================================================================\n";
  pout << "Current NPDM sweep position = "<< sweepPos+1 << " of " << endPos+1 << "\n";

  // Loop over NPDM operator patterns
  Npdm_patterns npdm_patterns( npdm_order_, sweepPos, endPos );

  if(wavefunctions.size()==2){
    Npdm_expectations npdm_expectations( spin_adaptation_, npdm_patterns, npdm_order_, wavefunctions.at(0), wavefunctions.at(1), big );
  //    for(auto pattern=npdm_patterns.ldr_cd_begin();pattern!=npdm_patterns.ldr_cd_end();++pattern){
  //      npdm_patterns.print_cd_string(pattern->at('l'));
  //      npdm_patterns.print_cd_string(pattern->at('d'));
  //      npdm_patterns.print_cd_string(pattern->at('r'));
  //      pout << endl;
  //    }
    loop_over_operator_patterns( npdm_patterns, npdm_expectations, big );
  }
  else{
    Npdm_expectations npdm_expectations( spin_adaptation_, npdm_patterns, npdm_order_, wavefunctions.at(0), wavefunctions.at(0), big );
  //    for(auto pattern=npdm_patterns.ldr_cd_begin();pattern!=npdm_patterns.ldr_cd_end();++pattern){
  //      npdm_patterns.print_cd_string(pattern->at('l'));
  //      npdm_patterns.print_cd_string(pattern->at('d'));
  //      npdm_patterns.print_cd_string(pattern->at('r'));
  //      pout << endl;
  //    }
    loop_over_operator_patterns( npdm_patterns, npdm_expectations, big );
  }

  // Print outs
#ifndef SERIAL
  if (mpigetrank() == 0) {
    int sum;
    reduce(world, DEBUG_CALL_GET_EXPECT[mpigetrank()], sum, std::plus<int>(), 0);
    pout << "NPDM calls to expectation engine " << sum << endl;
  } else {
    reduce(world, DEBUG_CALL_GET_EXPECT[mpigetrank()], std::plus<int>(), 0);
  }
  if (mpigetrank() == 0) {
    double sum;
    reduce(world, DEBUG_COMM_TIME[mpigetrank()], sum, std::plus<double>(), 0);
    pout << "NPDM mpi communications time " << sum << endl;
  } else {
    reduce(world, DEBUG_COMM_TIME[mpigetrank()], std::plus<double>(), 0);
  }
#endif
  pout << "NPDM compute elements time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
  pout << "===========================================================================================\n";
  
}

//===========================================================================================================================================================

}
}
