/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_DRIVER_HEADER_H
#define NPDM_DRIVER_HEADER_H

//FIXME remove includes?
#include <vector>
#include <multiarray.h>
#include "spinblock.h"
#include "wavefunction.h"
#include "BaseOperator.h"
#include "npdm_expectations.h"
#include "npdm_operator_wrappers.h"
#include "twopdm_container.h"
#include "threepdm_container.h"
#include "fourpdm_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Npdm_driver_base {
  public:
    Npdm_driver_base() {}
    virtual ~Npdm_driver_base() {}
    virtual void save_data() = 0;
    virtual void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) = 0;
};

//===========================================================================================================================================================

class Npdm_driver : public Npdm_driver_base {

  public:
    Npdm_driver(int order, Npdm_container& container) : npdm_order_(order), container_(container) {}
    ~Npdm_driver() {}
    void save_data() { container_.save_npdms(0,0); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos );

  private:
    int npdm_order_;
    Npdm_container& container_;

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

};
  
//===========================================================================================================================================================

class Twopdm_driver : public Npdm_driver_base {
  public:
    Twopdm_driver( int sites ) : container( Twopdm_container(sites) ), driver( Npdm_driver(2, container) ) {}
    Twopdm_container container;
    Npdm_driver driver;
    void save_data() { driver.save_data(); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
};

//===========================================================================================================================================================

class Threepdm_driver : public Npdm_driver_base {
  public:
    Threepdm_driver( int sites ) : container( Threepdm_container(sites) ), driver( Npdm_driver(3, container) ) {}
    Threepdm_container container;
    Npdm_driver driver;
    void save_data() { driver.save_data(); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
};

//===========================================================================================================================================================

class Fourpdm_driver : public Npdm_driver_base {
  public:
    Fourpdm_driver( int sites ) : container( Fourpdm_container(sites) ), driver( Npdm_driver(4, container) ) {}
    Fourpdm_container container;
    Npdm_driver driver;
    void save_data() { driver.save_data(); }
    void compute_npdm_elements( std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos ) 
      { driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos ); }
};

//===========================================================================================================================================================

}

#endif

