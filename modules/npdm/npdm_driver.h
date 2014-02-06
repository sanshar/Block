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
//#include "npdm_sparse_array.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Npdm_driver {

  public:
    explicit Npdm_driver(int order);
    virtual ~Npdm_driver() {};

    virtual void save_npdms(const int &i, const int &j) = 0;
    virtual void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos );

  protected:
    bool store_full_spin_array_;
    bool store_full_spatial_array_;

  private:
    int npdm_order_;

    void loop_over_operator_patterns( Npdm::Npdm_patterns& patterns, Npdm::Npdm_expectations& expectations, const SpinBlock& big );

    void do_inner_loop( const char inner, Npdm::Npdm_expectations & npdm_expectations, 
                        NpdmSpinOps_base & outerOps, NpdmSpinOps & innerOps, NpdmSpinOps & dotOps );

    void do_parallel_lhs_loop( const char inner, Npdm::Npdm_expectations & npdm_expectations,
                               NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps, bool skip );

    void loop_over_block_operators( const char inner, Npdm::Npdm_expectations & npdm_expectations, 
                                    NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, NpdmSpinOps & dotOps, bool lhsdot );

    int get_mpi_max_size( int my_size );
    bool broadcast_lhs( int lhs_size, int rhs_size );
    bool skip_this_mpi_rank( NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps );
    bool skip_parallel( NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps, bool lhsrhsdot );

    virtual void clear_sparse_arrays() = 0;
    virtual void update_full_spin_array() = 0;
    virtual void update_full_spatial_array() = 0;
    virtual void store_npdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements) = 0;

};
  
//===========================================================================================================================================================

}

#endif

