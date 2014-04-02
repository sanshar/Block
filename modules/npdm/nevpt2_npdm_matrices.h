/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NEVPT2_NPDM_HEADER
#define NEVPT2_NPDM_HEADER

#include "global.h"
#include "multiarray.h"

namespace SpinAdapted {
namespace Npdm {

//===========================================================================================================================================================

class Nevpt2_npdm {

  public:
    array_6d<double>  compute_EEE_matrix(array_2d<double>& onepdm, array_4d<double>& twopdm, array_6d<double>& threepdm );
    array_8d<double> compute_EEEE_matrix(array_2d<double>& onepdm, array_4d<double>& twopdm, array_6d<double>& threepdm, array_8d<double>& fourpdm );
    void compute_A16_matrix( array_8d<double>& eeee );
    void compute_A22_matrix( array_6d<double>& eee, array_8d<double>& eeee );

};

//===========================================================================================================================================================

}
}

#endif
