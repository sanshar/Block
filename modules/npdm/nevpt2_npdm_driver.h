/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_NPDM_DRIVER_H
#define NEVPT2_NPDM_DRIVER_H

//#define DEBUG_NEVPT2NPDM

#include "npdm_driver.h"
#include "nevpt2_container.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class Nevpt2_npdm_driver : public Npdm_driver_base {

  public:
    explicit Nevpt2_npdm_driver( int sites );
    ~Nevpt2_npdm_driver() {};
  
    void save_data( const int i, const int j );
    void compute_npdm_elements( std::vector<Wavefunction>& wavefunctions, const SpinBlock& big, int sweepPos, int endPos );
    void clear();

  private:
#ifdef DEBUG_NEVPT2NPDM
    // Build NPDMs first, construct A-matrices later
    Onepdm_driver onepdm_driver;
    Twopdm_driver twopdm_driver;
    Threepdm_driver threepdm_driver;
    Fourpdm_driver fourpdm_driver;
#else
    // Build A-matrices on the fly
    Nevpt2_container nevpt2_container;
    Npdm_driver onepdm_driver;
    Npdm_driver twopdm_driver;
    Npdm_driver threepdm_driver;
    Npdm_driver fourpdm_driver;
#endif
    void compute_matrices( const int i, const int j );

};

//===========================================================================================================================================================

}
}

#endif

