/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_HEADER
#define NEVPT2_HEADER

#include "npdm.h"

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void compute_EEEE_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& threepdm_driver, Fourpdm_driver& fourpdm_driver );
  void compute_A16_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& threepdm_driver, Fourpdm_driver& fourpdm_driver );

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}

#endif
