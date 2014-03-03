/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_NPDM_DRIVER_H
#define NEVPT2_NPDM_DRIVER_H

#include "npdm_driver.h"
#include "nevpt2_A16_container.h"

namespace SpinAdapted{

//===========================================================================================================================================================

class Nevpt2_npdm_driver : public Npdm_driver_base {

  public:
    Nevpt2_npdm_driver( int sites );
    ~Nevpt2_npdm_driver() {};
  
    void save_data();
    void compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos);

  private:

    Nevpt2_A16_matrix nevpt2_A16_matrix;
    Npdm_driver twopdm_driver;
    Npdm_driver threepdm_driver;
    Npdm_driver fourpdm_driver;

    void compute_matrices();
};

//===========================================================================================================================================================

}

#endif

