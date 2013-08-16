/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_DRIVER_HEADER_H
#define NPDM_DRIVER_HEADER_H

#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Npdm_driver {

public:
  Npdm_driver(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big);
  void compute_npdm_sweep(int state, int sweepPos, int endPos);

private:
  void npdm_loop_over_block_operators( std::vector<Npdm::CD> & lhs_cd,
                                       std::vector<Npdm::CD> & dot_cd,
                                       std::vector<Npdm::CD> & rhs_cd,
                                       array_4d<double> & npdm );

  std::vector< std::pair<bool, boost::shared_ptr<NpdmSpinOps>> > 
  get_all_mpi_ops( const bool local_skip, boost::shared_ptr<NpdmSpinOps>& local_ops, std::vector< boost::mpi::request > & reqs );

  void do_npdm_inner_loop( array_4d<double> & npdm,
                           boost::shared_ptr<NpdmSpinOps> & lhsOps,
                           boost::shared_ptr<NpdmSpinOps> & rhsOps,
                           boost::shared_ptr<NpdmSpinOps> & dotOps,
                           Npdm::Npdm_expectations & npdm_expectations );

  std::vector<Wavefunction> & wavefunctions_;
  const SpinBlock & big_; 

};

//===========================================================================================================================================================

}

#endif

