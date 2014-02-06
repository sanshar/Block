/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_PDM_DRIVER_H
#define NEVPT2_PDM_DRIVER_H

#include "npdm_driver.h"
#include "twopdm_container.h"
#include "threepdm_container.h"
#include "fourpdm_container.h"
//#include "contracted_fourpdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Nevpt2_pdm_driver : public Npdm_driver_base {

  public:
    Nevpt2_pdm_driver( int sites );
    ~Nevpt2_pdm_driver() {};
  
    void compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos);

  private:
    boost::shared_ptr<Twopdm_container> twopdm_container;
    boost::shared_ptr<Threepdm_container> threepdm_container;
    boost::shared_ptr<Fourpdm_container> fourpdm_container;

//    void update_A16_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& Threepdm_driver, Fourpdm_driver& fourpdm_driver );
//
//    // NEVPT2 driver only uses elements build by NPDM drivers, so doesn't need this function.
//    void store_npdm_elements(std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements) { assert(false); }
//    void clear_sparse_arrays() {}
//    void update_full_spin_array() {}
//    void update_full_spatial_array() {}

};

//===========================================================================================================================================================

}

#endif

