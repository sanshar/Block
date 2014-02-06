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

  array_6d<double> compute_EEE_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& threepdm_driver );
  array_8d<double> compute_EEEE_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& threepdm_driver, Fourpdm_driver& fourpdm_driver );
  void compute_A16_matrix( int dim, array_8d<double>& eeee );
  void compute_A22_matrix( int dim, array_6d<double>& eee, array_8d<double>& eeee );

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}

#endif
