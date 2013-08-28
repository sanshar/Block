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
// Base class
//===========================================================================================================================================================

class Npdm_driver {

  public:
  
    Npdm_driver() : npdm_order_(-1) { }
    Npdm_driver(int order) : npdm_order_(order) { }

    void compute_npdm_sweep(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int state, int sweepPos, int endPos);
    virtual void save_npdm_text(const int &i, const int &j) = 0;
    virtual void save_npdm_binary(const int &i, const int &j) = 0;
    virtual void save_spatial_npdm_text(const int &i, const int &j) = 0;
    virtual void save_spatial_npdm_binary(const int &i, const int &j) = 0;
    virtual void load_npdm_binary(const int &i, const int &j) = 0;
    virtual void resize_array(int dim) = 0;
    virtual void clear_array() = 0;
  
  protected:
  
    int npdm_order_;
  
    int get_mpi_max_lhs_size( int my_size );
    void do_npdm_inner_loop( Npdm::Npdm_expectations & npdm_expectations, NpdmSpinOps_base & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps );
    void npdm_loop_over_block_operators( Npdm::Npdm_expectations & npdm_expectations, NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps );
    std::vector<NpdmSpinOps_base> get_all_mpi_ops(const bool local_skip, NpdmSpinOps & local_ops, std::vector< boost::mpi::request > & reqs);
    virtual void assign_npdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements) = 0;
    virtual void accumulate_npdm() = 0;

};
  
//===========================================================================================================================================================

}

#endif

