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
#include "npdm_sparse_array.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Npdm_driver {

  public:
// FIXME flag to use full npdm array is hard-coded here; make input keyword
    Npdm_driver(bool full = true) : npdm_order_(-1), use_full_array_(full) {};
    Npdm_driver(int order, bool full = true) : npdm_order_(order), use_full_array_(full) {};
    virtual ~Npdm_driver() {};

    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos );
    void save_array(int i, int j);

  protected:
    bool use_full_array_;

  private:
    int npdm_order_;

    // We only store the nonredundant npdm matrix elements in these containers
//FIXME
//    Npdm_sparse_array sparse_spinorb_npdm_;
//    Npdm_sparse_array sparse_spatial_npdm_;
    Npdm_sparse_array sparse_array_;

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

    void save_full_array(int i, int j);
    void save_sparse_array(int i, int j);

    virtual void store_sparse_npdm_elements( std::vector< std::pair< std::vector<int>, double > > & elements) { sparse_array_.insert(elements); }

    virtual void store_npdm_elements( std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements) = 0;
    virtual void accumulate_npdm() = 0;

    virtual void save_npdm_text(const int &i, const int &j) = 0;
    virtual void save_npdm_binary(const int &i, const int &j) = 0;

    virtual void build_spatial_npdm(const int& i, const int& j) = 0;
    virtual void save_spatial_npdm_text(const int &i, const int &j) = 0;
    virtual void save_spatial_npdm_binary(const int &i, const int &j) = 0;

    virtual void load_npdm_binary(const int &i, const int &j) = 0;

};
  
//===========================================================================================================================================================

}

#endif

